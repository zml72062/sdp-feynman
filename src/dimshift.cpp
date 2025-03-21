#include "config.hpp"
#include "utils.hpp"
#include <fstream>

GiNaC::matrix config_parser::get_shift_to_upper_dim() {
    // map master integrals at dimension d-2 to linear combinations of
    // integrals at dimension d
    std::map<std::string, std::map<std::string, GiNaC::ex>> collection;
    std::map<std::string, int> master_to_id;
    int id = 0;
    GiNaC::symtab expansion;
    for (auto& master: master_table) {
        master_to_id[master] = id++;
        collection[master] = std::map<std::string, GiNaC::ex>();
        auto indices = split(master.c_str());
        auto termp = polynomial_iterator(symanzik_U), end = termp.end();
        for (; termp != end; ++termp) {
            std::vector<int> new_indices(feynman_params.nops(), 0);
            GiNaC::ex term = *termp;
            GiNaC::ex factor = 1;
            for (auto& i: effective_feynman_params) {
                int d = term.degree(feynman_params[i]);
                new_indices[i] = indices[i] + d;
                term = term.lcoeff(feynman_params[i]);
                for (int j = 0; j < d; j++) {
                    factor *= (indices[i] + j);
                }
            }
            factor *= term;
            auto repr = combine(new_indices);
            collection[master][repr] = factor;
            expansion[repr] = 0;
        }
    }

    // reduce integrals at dimension d to master integrals at dimension d
    std::ifstream ibp_result_file(ibp_result_filename);
    std::string ibp, current_key;
    std::map<std::string, std::map<std::string, GiNaC::ex>> storage;
    std::size_t _asterisk;
    bool exist = false;
    while (true) {
        ibp_result_file >> ibp;
        if (ibp_result_file.eof())
            break;
        if (ibp.find(integral_family) != std::string::npos) {
            if ((_asterisk = ibp.find('*')) == std::string::npos) { // an IBP head
                current_key = int_to_id(ibp);
                exist = (expansion.find(current_key) != expansion.end());
                if (exist)
                    storage[current_key] = std::map<std::string, GiNaC::ex>();
            } else if (exist) { // an IBP body
                auto current_integral = int_to_id(ibp);
                if (read_cache_exists(current_key, current_integral)) {
                    storage[current_key][current_integral] 
                        = read_ibp_simple(current_key, current_integral);
                } else {
                    if (working_subprocesses == max_subprocesses) {
                        read_subprocess_yield(false, [this, &storage](auto k, auto i) {
                            storage[k][i] = read_ibp_simple(k, i);
                        });
                    }
                    read_subprocess_work(current_key, current_integral, ibp.substr(_asterisk + 1));
                }
            }
        }
    }
    // add IBP entries for master integrals themselves
    for (auto& master: master_table) {
        storage[master] = std::map<std::string, GiNaC::ex>();
        storage[master][master] = 1;
    }
    while (working_subprocesses != 0)
        read_subprocess_yield(true, [this, &storage](auto k, auto i) {
            storage[k][i] = read_ibp_simple(k, i);
        });
    ibp_result_file.close();

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

