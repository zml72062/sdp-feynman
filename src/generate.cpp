#include "generate.hpp"
#include "utils.hpp"

static GiNaC::matrix get_quadratic_form(const GiNaC::ex& prefactor,
                                        const GiNaC::lst& terms) {
    int n = terms.nops();

    GiNaC::matrix out(n, n);
    int i = 0;
    for (auto p1 = terms.begin(); p1 != terms.end(); ++p1, i++) {
        int j = 0;
        for (auto p2 = terms.begin(); p2 != terms.end(); ++p2, j++) {
            out(i, j) = prefactor * (*p1) * (*p2);
        }
    }

    return out;
}

GiNaC::matrix polynomial_generator::generate(
    const std::string& prefactor, const std::vector<std::string>& terms) {
    
    GiNaC::ex prefactor_ex = (*parser)(prefactor);
    GiNaC::lst list_terms;
    for (auto& term: terms) {
        list_terms.append((*parser)(term));
    }

    return get_quadratic_form(prefactor_ex, list_terms);
}

GiNaC::matrix polynomial_generator::generate(
    const std::string& prefactor, int min_x_degree, int max_x_degree, int max_log_degree) {
    
    GiNaC::ex prefactor_ex = (*parser)(prefactor);
    GiNaC::ex x_factor = 1, log_factor = 1 + symbol_table["L"];
    GiNaC::lst rules, terms;
    GiNaC::symbol temp("temp");
    for (auto& key_value: symbol_table) {
        rules.append(key_value.second == 1);
        if (key_value.first == "L")
            continue;
        x_factor += (key_value.second * temp);
    }
    GiNaC::ex generating_poly = GiNaC::pow(x_factor, max_x_degree) 
                              * GiNaC::pow(log_factor, max_log_degree);
    if (min_x_degree > 0)
        generating_poly = generating_poly.diff(temp, min_x_degree);
    generating_poly = generating_poly.subs(temp == 1, GiNaC::subs_options::algebraic)
                                     .expand();

    auto iter = polynomial_iterator(generating_poly), end = iter.end();
    for (; iter != end; ++iter) {
        GiNaC::ex term = *iter;
        terms.append(term / term.subs(rules, GiNaC::subs_options::algebraic));
    }

    return get_quadratic_form(prefactor_ex, terms);
}

std::vector<GiNaC::matrix> polynomial_generator::generate_from_config() {
    auto ansatze = (*configp)["ansatze"].as<std::vector<YAML::Node>>();
    std::vector<GiNaC::matrix> results;
    for (auto& ansatz: ansatze) {
        auto parts = ansatz.as<std::vector<YAML::Node>>();
        if (parts.size() == 4) {
            // [prefactor, min_x_degree, max_x_degree, max_log_degree]
            auto prefactor = parts[0].as<std::string>();
            auto min_x_degree = parts[1].as<int>();
            auto max_x_degree = parts[2].as<int>();
            auto max_log_degree = parts[3].as<int>();
            results.push_back(generate(prefactor, min_x_degree, max_x_degree, max_log_degree));
        } else if (parts.size() == 2) {
            // [prefactor, [terms]]
            auto prefactor = parts[0].as<std::string>();
            auto terms = parts[1].as<std::vector<std::string>>();
            results.push_back(generate(prefactor, terms));
        }
    }

    return results;
}

