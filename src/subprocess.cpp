#include "config.hpp"
#include "utils.hpp"
#include <sys/wait.h>

void config_parser::subprocess_work(const std::string& key, const std::string& integral, const std::string& coefficient) {
    pid_t pid = fork();
    if (pid == 0) { // subprocess
        GiNaC::parser coefficient_reader(symbol_table);
        auto ex = coefficient_reader(coefficient);
        if (kinematics_numerics.nops() > 0) {
            ex = ex.subs(kinematics_numerics, GiNaC::subs_options::algebraic);
        }
        save_to_cache(key, integral, ex);
        exit(0);
    } else {
        working_subprocesses++;
        subprocess_map[pid] = std::make_pair(key, integral);
    }
}

void config_parser::mainprocess_work(const std::string& key, const std::string& integral) {
    auto current_integral = get(integral_table, integral, "I[", "]")
        .subs(master_values, GiNaC::subs_options::algebraic)
        .subs(kinematics_numerics, GiNaC::subs_options::algebraic);
    auto prefactor = get_prefactor(key, t, num_internals, symbol_table["d"],
                                   sector_designate, top_level_sector).second;
    auto coefficient = load_from_cache(key, integral);
    ibp_table[key] += (current_integral * coefficient * prefactor);
}

void config_parser::subprocess_yield(bool always_wait) {
    pid_t pid;
    int flag = always_wait ? 0 :
        ((working_subprocesses == max_subprocesses) ? 0 : WNOHANG);
    while ((pid = waitpid(-1, 0, flag)) > 0) {
        working_subprocesses--;
        auto key_integral = subprocess_map[pid];
        auto key = key_integral.first, integral = key_integral.second;
        mainprocess_work(key, integral);
        subprocess_map.erase(pid);
        flag = always_wait ? 0 :
            ((working_subprocesses == max_subprocesses) ? 0 : WNOHANG);
    }
}

