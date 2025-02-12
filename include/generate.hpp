#ifndef GENERATE_HPP
#define GENERATE_HPP

#include <ginac/ginac.h>
#include <yaml-cpp/yaml.h>

class polynomial_generator {
public:
    polynomial_generator(const GiNaC::lst& feynman_params,
                         const std::vector<int>& effective_feynman_params,
                         const GiNaC::symbol& log,
                         const YAML::Node& config)
        : configp(&config) {
        
        symbol_table["L"] = log;
        for (auto& param: effective_feynman_params) {
            auto var = feynman_params[param];
            symbol_table[GiNaC::ex_to<GiNaC::symbol>(var).get_name()] = var;
        }
        parser = new GiNaC::parser(symbol_table);
    }

    ~polynomial_generator() {
        delete parser;
    }

    /**
     * Generate a quadratic form of \{c_k\} of the following content:
     * 
     *      prefactor * (c_0 terms[0] + c_1 terms[1] + ...)^2
     * 
     * where `prefactor` and every `terms[k]` are polynomials of Feynman
     * parameters (\{x_i\}) and L = log(U^{L+1}/F^L).
     * 
     * @param prefactor a prefactor
     * @param terms a list of terms to be fully squared
     * 
     * @returns a symmetric square matrix representing the quadratic form
     */
    GiNaC::matrix generate(const std::string& prefactor, const std::vector<std::string>& terms);
    /**
     * Generate a quadratic form of \{c_k\} of the following content:
     * 
     *      prefactor * (c_0 terms[0] + c_1 terms[1] + ...)^2
     * 
     * where `prefactor` is a polynomial of Feynman parameters (\{x_i\}) and 
     * L = log(U^{L+1}/F^L). The `terms` enumerate all possible monomials
     * of degree at least `min_x_degree` and at most `max_x_degree` in the 
     * Feynman parameters (\{x_i\}), and of degree at most `max_log_degree` 
     * in L = log(U^{L+1}/F^L).
     * 
     * @param prefactor a prefactor
     * @param min_x_degree min degree in the Feynman parameters
     * @param max_x_degree max degree in the Feynman parameters
     * @param max_log_degree max degree in L = log(U^{L+1}/F^L)
     * 
     * @returns a symmetric square matrix representing the quadratic form
     */
    GiNaC::matrix generate(const std::string& prefactor, int min_x_degree, int max_x_degree, int max_log_degree);
    
    std::vector<GiNaC::matrix> generate_from_config();
private:
    GiNaC::symtab symbol_table;
    GiNaC::parser* parser;
    const YAML::Node* configp;
};


#endif
