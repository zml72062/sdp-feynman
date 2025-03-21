#include <fstream>
#include "utils.hpp"
#include "config.hpp"

void config_parser::read_ibps() {
    count_ibps();
    std::ifstream ibp_result_file(ibp_result_filename);
    std::string ibp, current_key;
    std::size_t _asterisk, counter = 0;
    bool fail = false;
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
                fail = !get_prefactor(current_key, t, num_internals, symbol_table["d"],
                                      sector_designate, top_level_sector).first;
                if (fail)
                    continue;
                get(integral_table, current_key, "I[", "]");
                ibp_table[current_key] = 0;
                std::cerr << "Processing the " << counter << "-th / " << ibp_count << " IBP relation" << "\r";
            } else if (!fail) { // an IBP body
                auto current_integral = int_to_id(ibp);
                if (read_cache_exists(current_key, current_integral)) {
                    read_mainprocess_work(current_key, current_integral);
                } else {
                    if (working_subprocesses == max_subprocesses) {
                        read_subprocess_yield(false, [this](auto k, auto i) { read_mainprocess_work(k, i); });
                    }
                    read_subprocess_work(current_key, current_integral, ibp.substr(_asterisk + 1));
                }
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
    while (working_subprocesses != 0)
        read_subprocess_yield(true, [this](auto k, auto i) { read_mainprocess_work(k, i); });
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
    int num_effective_ibps = ibp_table.size(), counter = 0;
    for (auto& key_value: ibp_table) {
        std::cerr << "Processing the " << ++counter << "-th / " << num_effective_ibps << " IBP relation" << "\r";
        if (expand_cache_exists(key_value.first)) {
            expand_mainprocess_work(key_value.first);
        } else {
            if (working_subprocesses == max_subprocesses) {
                expand_subprocess_yield(false);
            }
            expand_subprocess_work(key_value.first, key_value.second, rules, order, eps);
        }
    }
    while (working_subprocesses != 0)
        expand_subprocess_yield(true);
    END_TIME(expand_ibp);
    
    std::cerr << std::endl << "Done!" << std::endl;
    PRINT_TIME(expand_ibp);
}

void config_parser::expand_ibps() {
    expand_ibps(eps_order);
}

void config_parser::count_ibps() {
    std::ifstream ibp_result_file(ibp_result_filename);
    std::string ibp;
    int counter = 0;
    while (true) {
        ibp_result_file >> ibp;
        if (ibp_result_file.eof())
            break;
        if (ibp.find(integral_family) != std::string::npos
         && ibp.find('*') == std::string::npos)
            counter++;
    }
    ibp_result_file.close();
    ibp_count = counter;
}

void config_parser::dump_raw_ibps(std::ostream& out) {
    for (auto& ibp: ibp_table) {
        out << integral_table[ibp.first] << " = " << ibp.second << std::endl;
    }
}

void config_parser::dump_expanded_ibps(std::ostream& out) {
    int order_ = numeric_integral_table.size();
    for (int i = 0; i < order_; i++) {
        for (auto& ibp: numeric_ibp_table[i]) {
            out << numeric_integral_table[i][ibp.first] << " = "
                << ibp.second << std::endl;
        }
    }
}

