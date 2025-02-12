#include "solver.hpp"
#include "sdpa.hpp"
#include "utils.hpp"

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

void master_solver::solve_from(const std::vector<GiNaC::matrix>& matrices) {
    START_TIME(solve);
    fail = false;
    int num_blocks = matrices.size();
    int num_integrals = variables_to_solve.nops();

    std::vector<std::vector<GiNaC::matrix>> coefficients;
    std::vector<GiNaC::matrix> bias;
    GiNaC::lst zero_rules;
    for (auto& integral: variables_to_solve) {
        zero_rules.append(integral == 0);
    }
    for (int i = 0; i < num_integrals; i++) {
        coefficients.push_back(std::vector<GiNaC::matrix>(num_blocks));
        for (int j = 0; j < num_blocks; j++) {
            coefficients.back()[j] = GiNaC::ex_to<GiNaC::matrix>(
                matrices[j].diff(GiNaC::ex_to<GiNaC::symbol>(variables_to_solve[i]))
                           .subs(zero_rules, GiNaC::subs_options::algebraic)
            );
            if (all_zero(coefficients.back()[j])) {
                std::cerr << "Warning: " << variables_to_solve[i] << " cannot be determined" << std::endl;
            }
        }
    }
    for (int j = 0; j < num_blocks; j++) {
        bias.push_back(GiNaC::ex_to<GiNaC::matrix>(
            ((GiNaC::ex)matrices[j])
                .subs(zero_rules, GiNaC::subs_options::algebraic)
        ));
    }

    // instantiate an SDPA solver
    sdpa_interface solve(coefficients, bias, *configp);
    solve.solve();
    END_TIME(solve);
    PRINT_TIME(solve);

    if (solve.get_fail()) {
        fail = true;
        std::cerr << "SDPA failed at initialization!" << std::endl;
        return;
    }
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
}
