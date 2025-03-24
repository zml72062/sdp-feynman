#ifndef ASY_HPP
#define ASY_HPP

#include <ginac/ginac.h>
#include "utils.hpp"

class asy {
public:
    asy(const GiNaC::lst& feynman_params, 
        const std::vector<int>& effective_feynman_params,
        const GiNaC::ex& U, const GiNaC::ex& F,
        const GiNaC::symbol& kinematic_symbol)
    : feynman_paramsp(&feynman_params),
      effective_feynman_paramsp(&effective_feynman_params),
      UF_product((U * F).expand()), xp(&kinematic_symbol) { }

    GiNaC::matrix export_to_python() {
        std::vector<std::vector<int>> result;
        auto termp = polynomial_iterator(UF_product), end = termp.end();
        for (; termp != end; ++termp) {
            result.push_back(std::vector<int>());
            auto term = *termp;
            for (auto& i: *effective_feynman_paramsp) {
                result.back().push_back(term.degree((*feynman_paramsp)[i]));
                term = term.lcoeff((*feynman_paramsp)[i]);
            }
            result.back().push_back(term.degree(*xp));
        }
        
        if (result.size() == 0)
            return GiNaC::matrix(0, 0);
        
        int row = result.size(), col = result[0].size();
        GiNaC::matrix matrix(row, col);
        for (int i = 0; i < row; i++)
            for (int j = 0; j < col; j++)
                matrix(i, j) = result[i][j];
        return matrix;
    }

private:
    const GiNaC::lst* feynman_paramsp;
    const std::vector<int>* effective_feynman_paramsp;
    GiNaC::ex UF_product;
    const GiNaC::symbol* xp;
};


#endif // ASY_HPP
