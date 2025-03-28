#include "config.hpp"
#include "utils.hpp"
#include <fstream>

GiNaC::matrix config_parser::get_shift_to_upper_dim() {
    // map master integrals at dimension d-2 to linear combinations of
    // integrals at dimension d
    std::map<std::string, GiNaC::symtab> collection;
    std::map<std::string, int> master_to_id;
    int id = 0;
    GiNaC::symtab expansion;
    for (auto& master: master_table) {
        master_to_id[master] = id++;
        collection[master] = poly_to_terms(master, symanzik_U, &expansion);
    }

    // reduce integrals at dimension d to master integrals at dimension d
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

    return coeff_matrix;
}

GiNaC::matrix config_parser::get_shift_to_lower_dim() {
    return get_shift_to_upper_dim().inverse();
}

