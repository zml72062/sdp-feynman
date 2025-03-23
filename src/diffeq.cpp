#include "config.hpp"
#include "utils.hpp"
#include <fstream>

GiNaC::matrix config_parser::get_differential_equations(const GiNaC::symbol& symbol) {
    // map derivatives of master integrals at dimension d to linear
    // combinations of integrals at dimension d+2
    std::map<std::string, GiNaC::symtab> collection;
    std::map<std::string, int> master_to_id;
    int id = 0;
    GiNaC::symtab expansion;
    auto neg_deriv_F = -symanzik_F.subs(kinematics_numerics, GiNaC::subs_options::algebraic)
                                  .diff(symbol).expand();
    for (auto& master: master_table) {
        master_to_id[master] = id++;
        collection[master] = poly_to_terms(master, neg_deriv_F, &expansion);
    }

    // reduce integrals at dimension d+2 to master integrals at dimension d+2
    auto storage = read_selected_ibps(expansion);

    // fill in the matrix
    GiNaC::matrix coeff_matrix(id, id);
    for (int i = 0; i < id; i++) {
        for (int j = 0; j < id; j++) {
            coeff_matrix(i, j) = 0;
        }
    }
    for (auto& row: collection) {
        int r = master_to_id[row.first];
        for (auto& col: row.second) {
            if (storage.find(col.first) == storage.end()) {
                std::cerr << "The integral I[" << col.first << "] is absent in IBP!" << std::endl;
                return GiNaC::matrix(0, 0);
            }
            auto& col_map = storage[col.first];
            for (auto& key_value: col_map) {
                int c = master_to_id[key_value.first];
                auto v = key_value.second * col.second;
                coeff_matrix(r, c) += v;
            }
        }
    }

    // reduce master integrals at dimension d+2 to master integrals at 
    // dimension d
    auto dimshift = get_shift_to_lower_dim();
    if (dimshift.rows() == 0)
        return GiNaC::matrix(0, 0);
    
    GiNaC::symbol tempd("tempd");
    return GiNaC::ex_to<GiNaC::matrix>(
        ((GiNaC::ex)(coeff_matrix.mul(dimshift)))
            .subs(symbol_table["d"] == tempd + 2, GiNaC::subs_options::algebraic)
            .subs(tempd == symbol_table["d"], GiNaC::subs_options::algebraic)
    );
}

