#ifndef PARSE_HPP
#define PARSE_HPP

#include <ginac/ginac.h>
#include "generate.hpp"

class polynomial_parser {
public:
    polynomial_parser(const std::vector<std::string>& effective_master_table,
                      std::vector<GiNaC::symtab>& numeric_ibp_table,
                      const std::vector<int>& effective_feynman_params,
                      const GiNaC::lst& feynman_params,
                      const YAML::Node& config)
        : effective_master_tablep(&effective_master_table),
          numeric_ibp_tablep(&numeric_ibp_table),
          effective_feynman_paramsp(&effective_feynman_params),
          feynman_paramsp(&feynman_params),
          configp(&config), L("L") { }

    /**
     * Parse a polynomial of Feynman parameters (\{x_i\}) and
     * log(U^{L+1}/F^L) into a linear combination of master integrals.
     * 
     * @param polynomial the input polynomial whose arguments are 
     * Feynman parameters (\{x_i\}) and L = log(U^{L+1}/F^L)
     * @param matrix_input boolean flag indicating whether the input is
     * a matrix of polynomials (true) or a single polynomial (false)
     * 
     * @returns a `bool` value indicating whether parsing succeeds, and
     * the parsing result (if it does succeed)
     */
    std::pair<bool, GiNaC::ex> parse(const GiNaC::ex& polynomial, bool matrix_input);
    
    // Export internal data to a polynomial generator.
    polynomial_generator get_polynomial_generator() {
        return polynomial_generator(*feynman_paramsp,
                                    *effective_feynman_paramsp, 
                                    L, *configp);
    }
private:
    // these five pointers are owned by someone else
    const std::vector<std::string>* effective_master_tablep;
    std::vector<GiNaC::symtab>* numeric_ibp_tablep;
    const std::vector<int>* effective_feynman_paramsp;
    const GiNaC::lst* feynman_paramsp;
    const YAML::Node* configp;

    // "L" represents log(U^{L+1}/F^L)
    GiNaC::symbol L;
};


#endif

