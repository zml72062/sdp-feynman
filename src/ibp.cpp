#include <fstream>
#include "utils.hpp"
#include "config.hpp"

static std::pair<bool, GiNaC::ex> get_prefactor(const std::string& id, int t, int L,
                                                const GiNaC::ex& d,
                                                bool sector_designate,
                                                std::size_t top_level_sector) {
    auto indices = split(id.c_str());
    int len = indices.size();
    int total_a = 0;
    GiNaC::ex prefactor = 1;
    bool fail = false;
    for (int i = 0; i < len; i++) {
        if (sector_designate && !(top_level_sector & (1 << i)))
            continue;
        if (indices[i] <= 0) {
            fail = true;
            break;
        }
        prefactor *= GiNaC::tgamma(indices[i]);
        total_a += indices[i];
    }
    if (total_a < t)
        fail = true;
    if (fail)
        return std::make_pair(false, 1);
    
    for (int i = t; i < total_a; i++) {
        prefactor /= (i - d * L / 2);
    }
    return std::make_pair(true, prefactor);
}

void config_parser::read_ibps() {
    std::ifstream ibp_result_file(ibp_result_filename);
    std::string ibp, current_key;
    std::size_t _asterisk, counter = 0;
    bool fail = false;
    GiNaC::ex pre_factor;
    GiNaC::parser coefficient_reader(symbol_table);
    START_TIME(read_ibp);
    while (true) {
        ibp_result_file >> ibp;
        if (ibp_result_file.eof())
            break;
        if (ibp.find(integral_family) != std::string::npos) {
            if ((_asterisk = ibp.find('*')) == std::string::npos) { // an IBP head
                fail = false;
                counter++;
                current_key = int_to_id(ibp);
                auto output = get_prefactor(current_key, t, num_internals, symbol_table["d"],
                                            sector_designate, top_level_sector);
                fail = !output.first;
                pre_factor = output.second;
                if (fail)
                    continue;
                get(integral_table, current_key, "I[", "]");
                ibp_table[current_key] = 0;
                std::cerr << "Processing the " << counter << "-th IBP relation" << "\r";
            } else if (!fail) { // an IBP body
                auto current_integral = get(integral_table, int_to_id(ibp), "I[", "]")
                    .subs(master_values, GiNaC::subs_options::algebraic)
                    .subs(kinematics_numerics, GiNaC::subs_options::algebraic);
                auto coefficient = coefficient_reader(ibp.substr(_asterisk + 1));
                if (kinematics_numerics.nops() > 0) {
                    coefficient = coefficient.subs(kinematics_numerics, GiNaC::subs_options::algebraic);
                }
                ibp_table[current_key] += (current_integral * coefficient * pre_factor);
            }
        }
    }
    // add IBP entries for master integrals themselves
    for (auto& master: master_table) {
        auto rhs_integral = get(integral_table, master, "I[", "]")
            .subs(master_values, GiNaC::subs_options::algebraic)
            .subs(kinematics_numerics, GiNaC::subs_options::algebraic);
        auto output = get_prefactor(master, t, num_internals, symbol_table["d"],
                                    sector_designate, top_level_sector);
        if (!output.first)
            continue;
        ibp_table[master] = rhs_integral * output.second;
    }
    END_TIME(read_ibp);

    std::cerr << std::endl << "Done!" << std::endl;
    PRINT_TIME(read_ibp);
    ibp_result_file.close();
}

void config_parser::expand_ibps(int order) {
    START_TIME(expand_ibp);
    // generate integrals at different order
    numeric_integral_table = std::vector<GiNaC::symtab>(order + 1);
    for (int i = 0; i <= order; i++) {
        for (auto& key_value: integral_table) {
            get(numeric_integral_table[i], key_value.first, "I[", "]_" + std::to_string(i));
        }
    }
    // generate expansion rules
    GiNaC::symbol eps("eps");
    GiNaC::lst rules;
    rules.append(symbol_table["d"] == d0 - 2 * eps);
    for (auto& key_value: integral_table) {
        GiNaC::ex rhs = 0;
        for (int i = 0; i <= order; i++) {
            rhs += (numeric_integral_table[i][key_value.first] * GiNaC::pow(eps, i));
        }
        rules.append(key_value.second == rhs);
    }
    // generate IBP equations at different order
    numeric_ibp_table = std::vector<GiNaC::symtab>(order + 1);
    for (auto& key_value: ibp_table) {
        auto full = key_value.second.subs(rules, GiNaC::subs_options::algebraic);
        GiNaC::lst derivatives;
        int fail = 0;
        for (int i = 0; i <= order; i++) {
            try {
                derivatives.append(full.subs(eps == 0, GiNaC::subs_options::algebraic));
            } catch (GiNaC::pole_error& err) {
                fail = 1;
                break;
            }
            if (i != order)
                full = full.diff(eps);
        }
        if (fail)
            continue;
        for (int i = 0; i <= order; i++) {
            numeric_ibp_table[i][key_value.first] = derivatives[i] / GiNaC::tgamma(i + 1);
        }
    }
    END_TIME(expand_ibp);
    PRINT_TIME(expand_ibp);
}

void config_parser::expand_ibps() {
    expand_ibps(eps_order);
}

void config_parser::print_raw_ibps() {
    for (auto& ibp: ibp_table) {
        std::cout << integral_table[ibp.first] << " = " << ibp.second << std::endl;
    }
}

void config_parser::print_expanded_ibps() {
    int order_ = numeric_integral_table.size();
    for (int i = 0; i < order_; i++) {
        for (auto& ibp: numeric_ibp_table[i]) {
            std::cout << numeric_integral_table[i][ibp.first] << " = "
                      << ibp.second << std::endl;
        }
    }
}

