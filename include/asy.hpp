#ifndef ASY_HPP
#define ASY_HPP

#include <ginac/ginac.h>

class asy {
public:
    asy(const GiNaC::lst& feynman_params, 
        const std::vector<int>& effective_feynman_params,
        const GiNaC::ex& U, const GiNaC::ex& F,
        const GiNaC::symbol& kinematic_symbol, int L_, int d0_, 
        const GiNaC::symbol& d)
    : feynman_paramsp(&feynman_params),
      effective_feynman_paramsp(&effective_feynman_params),
      U_poly(U), F_poly(F), L(L_), d0(d0_),
      xp(&kinematic_symbol), dp(&d) { }

    void prepare(char* python_path) {
        python_work(python_path);
        compute_asymptotic_polys();
    }

    void try_integrate_at_boundary(const std::string& integral);

private:
    const GiNaC::lst* feynman_paramsp;
    const std::vector<int>* effective_feynman_paramsp;
    GiNaC::ex U_poly;
    GiNaC::ex F_poly;
    int L;
    int d0;
    const GiNaC::symbol* xp;
    const GiNaC::symbol* dp;
    GiNaC::matrix scaling_vectors;
    GiNaC::lst asymptotic_Us;
    GiNaC::lst asymptotic_Fs;

    GiNaC::matrix get_term_orders(const GiNaC::ex& polynomial);
    GiNaC::matrix export_to_python();
    void python_work(char* python_path);
    void compute_asymptotic_polys();
};


#endif // ASY_HPP
