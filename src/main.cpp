#include "config.hpp"
#include "utils.hpp"
#include <fstream>
#include <filesystem>

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <config_file.yaml>" << std::endl;
        exit(0);
    }
    
    config_parser configure(argv[1]);
    std::cout << "U polynomial = " << configure.U() << std::endl;
    std::cout << "F polynomial = " << configure.F() << std::endl;
    
#ifndef NO_GSL
    if (configure.will_check_euclidean)
        if (!configure.check_euclidean())
            exit(0);
#endif // NO_GSL
    
    configure.read_ibps();
    configure.expand_ibps();
    if (configure.will_dump_raw_ibps) {
        std::cerr << "Dumping raw IBPs ..." << std::endl;
        std::filesystem::create_directory("logs");
        std::ofstream raw_out(std::filesystem::path("logs").append("raw_ibps"));
        configure.dump_raw_ibps(raw_out);
        raw_out.close();
    }
    if (configure.will_dump_expanded_ibps) {
        std::cerr << "Dumping expanded IBPs ..." << std::endl;
        std::filesystem::create_directory("logs");
        std::ofstream expanded_out(std::filesystem::path("logs").append("expanded_ibps"));
        configure.dump_expanded_ibps(expanded_out);
        expanded_out.close();
    }

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
