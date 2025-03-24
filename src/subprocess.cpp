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
    int order = lst.nops();
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

void config_parser::generate_subprocess_work(int integral, int block, time_t timestamp, const GiNaC::matrix& raw_matrix, const GiNaC::lst& rules, const GiNaC::symbol& integral_symbol) {
    pid_t pid = fork();
    if (pid == 0) { // subprocess
        if (integral != -1) // coefficient
            save_to_generate_cache(integral, block, timestamp, GiNaC::ex_to<GiNaC::matrix>(raw_matrix.diff(integral_symbol).subs(rules, GiNaC::subs_options::algebraic)));
        else // bias
            save_to_generate_cache(integral, block, timestamp, GiNaC::ex_to<GiNaC::matrix>(((GiNaC::ex)raw_matrix).subs(rules, GiNaC::subs_options::algebraic)));
        exit(0);
    } else {
        working_subprocesses++;
        generate_subprocess_map[pid] = std::make_pair(integral, block);
    }
}

void config_parser::generate_mainprocess_work(int integral, int block, time_t timestamp, std::vector<std::vector<GiNaC::matrix>>* coefficient, std::vector<GiNaC::matrix>* bias) {
    if (integral == -1) { // bias
        (*bias)[block] = load_from_generate_cache(integral, block, timestamp);
    } else { // coefficient
        (*coefficient)[integral][block] = load_from_generate_cache(integral, block, timestamp);
    }
}

void config_parser::generate_subprocess_yield(bool always_wait, time_t timestamp, std::vector<std::vector<GiNaC::matrix>>* coefficient, std::vector<GiNaC::matrix>* bias) {
    pid_t pid;
    int flag = always_wait ? 0 :
        ((working_subprocesses == max_subprocesses) ? 0 : WNOHANG);
    while ((pid = waitpid(-1, 0, flag)) > 0) {
        working_subprocesses--;
        auto pair = generate_subprocess_map[pid];
        int integral = pair.first, block = pair.second;
        generate_mainprocess_work(integral, block, timestamp, coefficient, bias);
        generate_subprocess_map.erase(pid);
        flag = always_wait ? 0 :
            ((working_subprocesses == max_subprocesses) ? 0 : WNOHANG);
    }
}

