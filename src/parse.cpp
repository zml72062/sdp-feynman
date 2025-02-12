#include "utils.hpp"
#include "parse.hpp"


std::pair<bool, GiNaC::ex> polynomial_parser::parse(const GiNaC::ex& polynomial, bool matrix_input) {
    if (matrix_input) {
        auto matrix = GiNaC::ex_to<GiNaC::matrix>(polynomial);
        int r = matrix.rows(), c = matrix.cols();
        GiNaC::matrix result(r, c);
        bool fail = false;
        for (int i = 0; i < r; i++) {
            if (fail)
                break;
            for (int j = 0; j < c; j++) {
                if (fail)
                    break;
                auto output = parse(matrix(i, j), false);
                if (!output.first)
                    fail = true;
                result(i, j) = output.second;
            }
        }
        if (fail)
            return std::make_pair(false, result);
        
        return std::make_pair(true, result);
    }
    GiNaC::ex expanded_polynomial = polynomial.expand();

    // parse each term of the expanded polynomial
    GiNaC::ex generated_polynomial = 0;
    auto termp = polynomial_iterator(expanded_polynomial), end = termp.end();
    try {
        for (; termp != end; ++termp) {
            std::vector<int> integral_indices(feynman_paramsp->nops(), 0);
            int log_power;
            GiNaC::ex term = *termp;
            for (auto& i: *effective_feynman_paramsp) {
                integral_indices[i] = term.degree((*feynman_paramsp)[i]) + 1;
                term = term.lcoeff((*feynman_paramsp)[i]);
            }
            std::string indices_repr = combine(integral_indices);

            log_power = term.degree(L);
            GiNaC::ex coeff = term.lcoeff(L) * GiNaC::tgamma(log_power + 1);

            generated_polynomial += (coeff * numeric_ibp_tablep->at(log_power).at(indices_repr));
        }
    } catch (std::out_of_range& error) { // thrown by at()
        std::cerr << "The term " << *termp << " corresponds to no entry "
                  << "in the IBP table!" << std::endl;
        std::cerr << "Polynomial parser exiting..." << std::endl;
        return std::make_pair(false, (GiNaC::ex)0);
    }

    return std::make_pair(true, generated_polynomial);
}

