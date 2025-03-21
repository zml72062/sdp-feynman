#include "config.hpp"
#include "utils.hpp"
#include <sys/wait.h>
#include <functional>

void config_parser::read_subprocess_work(const std::string& key, const std::string& integral, const std::string& coefficient) {
    pid_t pid = fork();
    if (pid == 0) { // subprocess
        GiNaC::parser coefficient_reader(symbol_table);
        auto ex = coefficient_reader(coefficient);
        if (kinematics_numerics.nops() > 0) {
            ex = ex.subs(kinematics_numerics, GiNaC::subs_options::algebraic);
        }
        save_to_read_cache(key, integral, ex);
        exit(0);
    } else {
        working_subprocesses++;
        read_subprocess_map[pid] = std::make_pair(key, integral);
    }
}

void config_parser::read_mainprocess_work(const std::string& key, const std::string& integral) {
    auto current_integral = get(integral_table, integral, "I[", "]")
        .subs(master_values, GiNaC::subs_options::algebraic)
        .subs(kinematics_numerics, GiNaC::subs_options::algebraic);
    auto prefactor = get_prefactor(key, t, num_internals, symbol_table["d"],
                                   sector_designate, top_level_sector).second;
    ibp_table[key] += (current_integral * prefactor * read_ibp_simple(key, integral));
}

void config_parser::read_subprocess_yield(bool always_wait, const std::function<void(const std::string&, const std::string&)>& callback) {
    pid_t pid;
    int flag = always_wait ? 0 :
        ((working_subprocesses == max_subprocesses) ? 0 : WNOHANG);
    while ((pid = waitpid(-1, 0, flag)) > 0) {
        working_subprocesses--;
        auto key_integral = read_subprocess_map[pid];
        auto key = key_integral.first, integral = key_integral.second;
        callback(key, integral);
        read_subprocess_map.erase(pid);
        flag = always_wait ? 0 :
            ((working_subprocesses == max_subprocesses) ? 0 : WNOHANG);
    }
}

void config_parser::expand_subprocess_work(const std::string& key, const GiNaC::ex& coefficient, const GiNaC::lst& rules, int order, const GiNaC::symbol& eps) {
    pid_t pid = fork();
    if (pid == 0) { // subprocess
        auto full = coefficient.subs(rules, GiNaC::subs_options::algebraic);
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
            save_to_expand_cache(key, GiNaC::lst());
        else
            save_to_expand_cache(key, derivatives);
        exit(0);
    } else {
        working_subprocesses++;
        expand_subprocess_map[pid] = key;
    }
}

void config_parser::expand_mainprocess_work(const std::string& key) {
    auto lst = GiNaC::ex_to<GiNaC::lst>(load_from_expand_cache(key));
    auto order = lst.nops();
    if (order == 0)
        return;
    
    for (int i = 0; i < order; i++) {
        numeric_ibp_table[i][key] = lst[i] / GiNaC::tgamma(i + 1);
    }
}

void config_parser::expand_subprocess_yield(bool always_wait) {
    pid_t pid;
    int flag = always_wait ? 0 :
        ((working_subprocesses == max_subprocesses) ? 0 : WNOHANG);
    while ((pid = waitpid(-1, 0, flag)) > 0) {
        working_subprocesses--;
        auto key = expand_subprocess_map[pid];
        expand_mainprocess_work(key);
        expand_subprocess_map.erase(pid);
        flag = always_wait ? 0 :
            ((working_subprocesses == max_subprocesses) ? 0 : WNOHANG);
    }
}

