#include "solver.hpp"
#include "sdpa.hpp"
#include "utils.hpp"
#include "config.hpp"
#include <time.h>
#include <fstream>
#include <filesystem>

static bool all_zero(const GiNaC::matrix& matrix) {
    int r = matrix.rows(), c = matrix.cols();
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            if (!(matrix(i, j) == 0))
                return false;
        }
    }
    return true;
}

void master_solver::solve_from(const std::vector<GiNaC::matrix>& matrices, void* config_parserp) {
    START_TIME(solve);
    fail = false;
    int num_blocks = matrices.size();
    int num_integrals = variables_to_solve.nops();
    config_parser* configurep = (config_parser*)config_parserp;
    time_t now;
    time(&now);
    std::cerr << "Current time: " << now << std::endl;

    std::filesystem::create_directory("logs");
    std::ofstream variables_out(std::filesystem::path("logs").append("variables_to_solve"));
    variables_out << variables_to_solve << std::endl;
    variables_out.close();

    std::vector<std::vector<GiNaC::matrix>> coefficients;
    for (int i = 0; i < num_integrals; i++) {
        coefficients.push_back(std::vector<GiNaC::matrix>(num_blocks));
    }
    std::vector<GiNaC::matrix> bias(num_blocks);
    GiNaC::lst zero_rules;
    for (auto& integral: variables_to_solve) {
        zero_rules.append(integral == 0);
    }

    std::cerr << "Start generating SDP problem ..." << std::endl;
    int total_matrices = (num_integrals + 1) * num_blocks, cnt = 0;
    for (int i = 0; i < num_integrals; i++) {
        for (int j = 0; j < num_blocks; j++) {
            if (configurep->working_subprocesses == configurep->max_subprocesses) {
                configurep->generate_subprocess_yield(false, now, &coefficients, &bias);
            }
            std::cerr << "Processing " << ++cnt << "-th / " << total_matrices << " matrix" << "\r";
            configurep->generate_subprocess_work(i, j, now, matrices[j], zero_rules, GiNaC::ex_to<GiNaC::symbol>(variables_to_solve[i]));
        }
    }
    for (int j = 0; j < num_blocks; j++) {
        if (configurep->working_subprocesses == configurep->max_subprocesses) {
            configurep->generate_subprocess_yield(false, now, &coefficients, &bias);
        }
        std::cerr << "Processing " << ++cnt << "-th / " << total_matrices << " matrix" << "\r";
        configurep->generate_subprocess_work(-1, j, now, matrices[j], zero_rules, GiNaC::symbol("#####placeholder#####"));
    }
    while (configurep->working_subprocesses != 0)
        configurep->generate_subprocess_yield(true, now, &coefficients, &bias);

    std::cerr << std::endl;

    for (int i = 0; i < num_integrals; i++) {
        bool zero = true;
        for (int j = 0; j < num_blocks; j++) {
            zero &= all_zero(coefficients[i][j]);
        }
        if (zero) {
            std::cerr << "Warning: " << variables_to_solve[i] << " cannot be determined" << std::endl;
        }
    }
    std::cerr << "Finished generating SDP problem!" << std::endl;

    // instantiate an SDPA solver
    sdpa_interface solve(coefficients, bias, *configp);
    if (solve.get_fail()) {
        fail = true;
        std::cerr << "SDPA failed at initialization!" << std::endl;
        return;
    }

    if (will_dump) {
        std::cerr << "Dumping symbolic SDP ..." << std::endl;
        std::filesystem::create_directory("logs");
        std::ofstream sdp_out(std::filesystem::path("logs").append("symbolic_sdp"));
        solve.dump(sdp_out);
        sdp_out.close();
    }
    solve.solve();
    END_TIME(solve);
    PRINT_TIME(solve);

#ifndef NO_SDPA_LIB
    auto& result = solve.get_result();
    if (result.phase != SDPA::PhaseType::pdOPT) {
        fail = true;
        std::cerr << "SDPA failed after trying to solve!" << std::endl;
        return;
    }

    auto& result_vec = result.x_vec;
    std::cout << "Maximized minimum eigenvalue: " << -result_vec.back() << std::endl;
    if (result_vec.back() > positivity_threshold) {
        fail = true;
        std::cerr << "SDPA asserts that positivity constraints are infeasible!" << std::endl;
        return;
    }

    computed_values = GiNaC::lst();
    for (int i = 0; i < num_integrals; i++) {
        computed_values.append(variables_to_solve[i] == result_vec[i]);
    }
    fail = false;
#endif // NO_SDPA_LIB
}
