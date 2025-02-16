#include <fstream>
#include <filesystem>
#include "utils.hpp"
#include "config.hpp"

config_parser::config_parser(const char* _config_file_name) {
    START_TIME(initialize);
    config_file = YAML::LoadFile(_config_file_name);
    integral_family = config_file["integralfamily"]["name"].as<std::string>();
    kira_dir = config_file["kiradir"].as<std::string>();
    ibp_result_filename = std::filesystem::path(kira_dir)
        .append("results")
        .append(integral_family)
        .append(config_file["kirafile"].as<std::string>());
    t = config_file["t"].as<int>();
    d0 = config_file["d0"].as<int>();
    eps_order = config_file["eps_order"].as<int>();
    
    read_internals();
    read_externals();
    read_invariants();
    read_substitution_rules();
    read_propagators();
    read_symbols();
    read_kinematics_numerics();
    read_masters();
    compute_symanzik();
    read_master_values();

    END_TIME(initialize);
    PRINT_TIME(initialize);
}

void config_parser::read_internals() {
    auto internals = config_file["integralfamily"]["internals"].as<std::vector<std::string>>();
    put_raw_symbols(internals);
    num_internals = internals.size();
}

void config_parser::read_externals() {
    put_raw_symbols(config_file["integralfamily"]["externals"].as<std::vector<std::string>>());
}

void config_parser::read_invariants() {
    auto invariants = config_file["integralfamily"]["kinematic_invariants"].as<std::vector<YAML::Node>>();
    for(auto& invariant: invariants) {
        get(symbol_table, invariant.as<std::vector<std::string>>()[0]);
    }
}

void config_parser::read_substitution_rules() {
    GiNaC::parser parser(symbol_table);

    auto scalarproduct_rules = config_file["integralfamily"]["scalarproduct_rules"].as<std::vector<YAML::Node>>();
    for (auto& rule: scalarproduct_rules) {
        auto rule_vector = rule.as<std::vector<YAML::Node>>();
        auto rule_head = rule_vector[0].as<std::vector<std::string>>();
        substitution_rules.append(parser(rule_head[0]) * parser(rule_head[1]) == parser(rule_vector[1].as<std::string>()));
    }
}

void config_parser::read_propagators() {
    GiNaC::parser parser(symbol_table);
    
    auto _propagators = config_file["integralfamily"]["propagators"].as<std::vector<YAML::Node>>();
    for (auto& propagator: _propagators) {
        auto propagator_content = propagator.as<std::vector<std::string>>();
        propagators.append(
            (GiNaC::pow(parser(propagator_content[0]), 2) 
                - parser(propagator_content[1]))
            .expand()
            .subs(substitution_rules, GiNaC::subs_options::algebraic)
        );
    }
}

void config_parser::read_symbols() {
    std::ifstream symbol_file;
    std::filesystem::path symbol_path = std::filesystem::path(kira_dir)
        .append("sectormappings").append("variables");
    symbol_file.open(symbol_path);
    std::string symbol;
    while (true) {
        symbol_file >> symbol;
        if (symbol_file.eof())
            break;
        get(symbol_table, symbol);
    }
    symbol_file.close();
}

void config_parser::read_kinematics_numerics() {
    if (has_non_null_key(config_file, "kinematics_numerics")) {
        auto kinematics_numerics_vector = config_file["kinematics_numerics"]
            .as<std::vector<YAML::Node>>();
        GiNaC::parser value_parser; // a parser for numeric values
        for (auto& v: kinematics_numerics_vector) {
            auto key_value = v.as<std::vector<std::string>>();
            kinematics_numerics.append(
                get(symbol_table, key_value[0])
                    == GiNaC::ex_to<GiNaC::numeric>(value_parser(key_value[1]))
            );
        }
    }
}

void config_parser::read_masters() {
    std::ifstream masters_file;
    std::filesystem::path masters_path = std::filesystem::path(kira_dir)
        .append("results")
        .append(integral_family)
        .append("masters.final");
    masters_file.open(masters_path);
    std::string master;
    while (true) {
        masters_file >> master;
        if (masters_file.eof())
            break;
        if (master.find(integral_family) != std::string::npos) {
            std::string id = int_to_id(master);
            master_table.push_back(id);
            get(integral_table, id, "I[", "]");
        }
    }
    masters_file.close();
}

void config_parser::compute_symanzik() {
    std::size_t num_propagators = propagators.nops();
    for (std::size_t i = 0; i < num_propagators; i++) {
        feynman_params.append(GiNaC::symbol("x" + std::to_string(i)));
    }
    
    GiNaC::ex denominator = 0;
    sector_designate = has_non_null_key(config_file["integralfamily"], "top_level_sector");
    top_level_sector = 0;
    if (sector_designate)
        top_level_sector = config_file["integralfamily"]["top_level_sector"].as<std::size_t>();
    for (std::size_t i = 0; i < num_propagators; i++) {
        if (sector_designate && !(top_level_sector & (1 << i)))
            continue;
        denominator += (feynman_params[i] * propagators[i]);
        effective_feynman_params.push_back(i);
    }
    GiNaC::lst all_zero;
    auto internals = config_file["integralfamily"]["internals"].as<std::vector<std::string>>();
    GiNaC::lst internal_symbols;
    for (auto& internal: internals) {
        auto momentum = get(symbol_table, internal);
        internal_symbols.append(momentum);
        all_zero.append(momentum == 0);
    }

    GiNaC::ex J = -denominator.subs(all_zero, GiNaC::subs_options::algebraic);
    GiNaC::matrix V(num_internals, 1);
    GiNaC::matrix M(num_internals, num_internals);
    for (std::size_t i = 0; i < num_internals; i++) {
        GiNaC::ex deriv = denominator.diff(GiNaC::ex_to<GiNaC::symbol>(internal_symbols[i]));
        V(i, 0) = -deriv.subs(all_zero, GiNaC::subs_options::algebraic) / 2;
        for (std::size_t j = 0; j < num_internals; j++) {
            M(i, j) = deriv.diff(GiNaC::ex_to<GiNaC::symbol>(internal_symbols[j])) / 2;
        }
    }

    symanzik_U = M.determinant().expand();
    symanzik_F = (V.transpose().mul(adjugate(M)).mul(V)(0, 0) + J * symanzik_U
                    ).expand().subs(substitution_rules, GiNaC::subs_options::algebraic).expand();
}

void config_parser::read_master_values() {
    std::set<std::string> known_masters;
    GiNaC::parser parser(symbol_table);
    auto _values = config_file["master_values"].as<std::vector<YAML::Node>>();
    for (auto& _value: _values) {
        auto key_value = _value.as<std::vector<YAML::Node>>();
        auto keys = key_value[0].as<std::vector<YAML::Node>>();
        int num_keys = keys.size();
        std::string repr;
        for (int i = 0; i < num_keys; i++) {
            repr += keys[i].as<std::string>();
            if (i != num_keys - 1)
                repr += ",";
        }
        known_masters.insert(repr);
        master_values.append(integral_table[repr] == parser(key_value[1].as<std::string>()));
    }
    for (auto& master: master_table) {
        if (known_masters.find(master) == known_masters.end())
            effective_master_table.push_back(master);
    }
}

bool config_parser::check_euclidean(unsigned long seed, int trials) {
    std::cerr << "Euclidean checks start..." << std::endl;
    START_TIME(check);
    auto numeric_F = symanzik_F.subs(kinematics_numerics, GiNaC::subs_options::algebraic);
    int n = effective_feynman_params.size();
    random_feynman_params rng(n, seed);
    std::vector<double> random_instance(n);
    min_log_value = max_log_threshold;
    max_log_value = min_log_threshold;
    for (int i = 0; i < trials; i++) {
        rng(random_instance);
        GiNaC::lst rules;
        for (int j = 0; j < n; j++)
            rules.append(feynman_params[effective_feynman_params[j]] == random_instance[j]);
        double evaluated_F = 0;
        try {
            evaluated_F = to_double(numeric_F.subs(rules, GiNaC::subs_options::algebraic));
        } catch (...) {
            std::cerr << "Incomplete numerics for kinematics!" << std::endl;
            return false;
        }
        if (evaluated_F <= 0) {
            std::cerr << std::endl << "Euclidean check failed!" << std::endl;
            return false;
        }
        double evaluated_U = to_double(symanzik_U.subs(rules, GiNaC::subs_options::algebraic));
        double value = std::pow(evaluated_U, num_internals + 1) / std::pow(evaluated_F, num_internals);
        if (value > max_log_value) {
            max_log_value = value;
            max_log_point = random_instance;
        }
        if (value < min_log_value) {
            min_log_value = value;
            min_log_point = random_instance;
        }
        std::cerr << i + 1 << " / " << trials << " trials..." << "\r";
    }
    END_TIME(check);
    std::cerr << std::endl << "Euclidean check succeeded!" << std::endl;
    PRINT_TIME(check);
    std::cout << "Detected max(U^{L+1}/F^L) = " << max_log_value << " at ";
    GiNaC::lst max_rule;
    for (int i = 0; i < n; i++) {
        max_rule.append(feynman_params[effective_feynman_params[i]] == max_log_point[i]);
    }
    std::cout << max_rule << std::endl;
    std::cout << "Detected min(U^{L+1}/F^L) = " << min_log_value << " at ";
    GiNaC::lst min_rule;
    for (int i = 0; i < n; i++) {
        min_rule.append(feynman_params[effective_feynman_params[i]] == min_log_point[i]);
    }
    std::cout << min_rule << std::endl;
    return true;
}

