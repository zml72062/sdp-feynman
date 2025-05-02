#include "config.hpp"
#include "utils.hpp"
#include <fstream>
#include <filesystem>

#define NORMAL(mat) ((GiNaC::ex)(mat)).normal()

int main(int argc, const char** argv) {
    if (argc == 1) {
        config_parser configure_diff("additional_examples/irexample_configs/diff.yaml");
        configure_diff.read_ibps();
        configure_diff.expand_ibps();
        std::cout << "Differential equation w.r.t. auxiliary mass parameter x is:" << std::endl;
        std::cout << "(d/dx) I(x) = A(x) I(x), where A(x) is" << std::endl;
        std::cout << NORMAL(configure_diff.get_differential_equations()) << std::endl;
        exit(0);
    }

    if (argc == 2) {
        std::string param(argv[1]);
        std::cout << "Now evaluate epsilon expansions at x=" << param << std::endl;
        
        GiNaC::Digits = 64;
        config_parser configure(("additional_examples/irexample_configs/eval" + param + ".yaml").c_str());
        configure.read_ibps();
        configure.expand_ibps();

        auto parser = configure.get_polynomial_parser();
        auto generator = parser.get_polynomial_generator();
        auto polynomials = generator.generate_from_config();

        START_TIME(parse_polynomials);
        std::vector<GiNaC::matrix> matrices;
        for (auto& polynomial: polynomials) {
            auto parser_output = parser.parse(polynomial, true);
            if (parser_output.first)
                matrices.push_back(GiNaC::ex_to<GiNaC::matrix>(parser_output.second));
        }
        END_TIME(parse_polynomials);
        PRINT_TIME(parse_polynomials);

        if (matrices.size() == 0) {
            std::cerr << "No available positivity constraints!" << std::endl;
            std::cerr << "Exiting..." << std::endl;
            exit(0);
        }

        auto solver = configure.get_solver();
        solver.solve_from(matrices, &configure);

#ifndef NO_SDPA_LIB
        if (solver.get_fail())
            exit(0);

        std::cout << "Computed master integral values are:" << std::endl;
        std::cout << solver.get_result() << std::endl;
#endif // NO_SDPA_LIB
    }
}
