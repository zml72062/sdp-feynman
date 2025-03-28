#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <yaml-cpp/yaml.h>
#include <ginac/ginac.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>
#include "parse.hpp"
#include "solver.hpp"
#include "utils.hpp"
#include "asy.hpp"


class config_parser {
public:
    /**
     * Initialize the parser from YAML file.
     * @param _config_file_name YAML config file name
     */
    config_parser(const char* _config_file_name);

#ifndef NO_GSL
    /**
     * Check whether the given kinematics lies in the Euclidean region.
     * Numerics for *all* kinematics should be given in YAML config file.
     * @param seed random seed for generating numeric Feynman parameters
     * @param trials number of trials
     */
    bool check_euclidean(unsigned long seed = 0, int trials = 65536);
#endif // NO_GSL

    // IBP related methods
    void read_ibps();
    void expand_ibps(int order);
    void expand_ibps();
    void dump_raw_ibps(std::ostream& out);
    void dump_expanded_ibps(std::ostream& out);
    std::map<std::string, GiNaC::symtab> read_selected_ibps(const GiNaC::symtab& integrals);

    // dimensional shifting relations and differential equations
    GiNaC::matrix get_shift_to_upper_dim();
    GiNaC::matrix get_shift_to_lower_dim();
    GiNaC::matrix get_differential_equations(const GiNaC::symbol& symbol);

    // Export internal data to a polynomial parser.
    polynomial_parser get_polynomial_parser() {
        return polynomial_parser(effective_master_table, 
                                 numeric_ibp_table,
                                 effective_feynman_params,
                                 feynman_params,
                                 config_file);
    }

    // Export internal data to a master solver.
    master_solver get_solver() {
        return master_solver(effective_master_table,
                             numeric_integral_table,
                             config_file, will_dump_symbolic_sdp);
    }

    // Export internal data to asy.
    asy get_asy(const GiNaC::symbol& kinematic_symbol) {
        return asy(feynman_params,
                   effective_feynman_params,
                   U(), F(), kinematic_symbol,
                   num_internals, d0, 
                   GiNaC::ex_to<GiNaC::symbol>(symbol_table["d"]));
    }

    const GiNaC::symbol* get_diff_variablep() {
        if (diff_variable.get_name().find("symbol") == 0)
            return nullptr;
        else
            return &diff_variable;
    }

    const std::string& family_name() {
        return integral_family;
    }

    const GiNaC::ex& U() {
        return symanzik_U;
    }

    const GiNaC::ex& F() {
        return symanzik_F;
    }

    // Options
    bool will_check_euclidean;
    bool will_dump_raw_ibps;
    bool will_dump_expanded_ibps;
    bool will_dump_symbolic_sdp;

    friend class master_solver;
private:
    YAML::Node config_file;
    // name of the family of integrals
    std::string integral_family;
    // directory where Kira places its output
    std::string kira_dir;
    // file name of Kira IBP reduction result
    std::string ibp_result_filename;
    
    GiNaC::symtab symbol_table;
    GiNaC::symtab integral_table;
    GiNaC::symtab ibp_table;
    std::vector<std::string> master_table;
    GiNaC::lst kinematics_numerics;
    std::vector<GiNaC::symtab> numeric_integral_table;
    std::vector<GiNaC::symtab> numeric_ibp_table;
    GiNaC::lst master_values;
    std::vector<std::string> effective_master_table;

    GiNaC::lst substitution_rules;
    GiNaC::lst propagators;

    GiNaC::lst feynman_params;
    std::vector<int> effective_feynman_params;
    GiNaC::ex symanzik_U;
    GiNaC::ex symanzik_F;

    std::vector<double> min_log_point;
    std::vector<double> max_log_point;
    const double min_log_threshold = 1e-9;
    const double max_log_threshold = 1e9;
    double min_log_value;
    double max_log_value;

    int ibp_count;
    int d0;
    int t;
    int eps_order;
    int num_internals;
    bool sector_designate;
    std::size_t top_level_sector;
    GiNaC::symbol diff_variable;

    void read_internals();
    void read_externals();
    void read_invariants();
    void read_substitution_rules();
    void read_propagators();
    void read_symbols();
    void read_kinematics_numerics();
    void read_masters();
    void compute_symanzik();
    void read_master_values();
    void count_ibps();

    // cache management
    std::string cache_dir;
    bool read_cache_exists(const std::string& key, const std::string& integral);
    GiNaC::ex load_from_read_cache(const std::string& key, const std::string& integral);
    void save_to_read_cache(const std::string& key, const std::string& integral, const GiNaC::ex& coefficient);
    GiNaC::ex read_ibp_simple(const std::string& key, const std::string& integral);
    bool expand_cache_exists(const std::string& key);
    GiNaC::ex load_from_expand_cache(const std::string& key);
    void save_to_expand_cache(const std::string& key, const GiNaC::ex& coefficient);
    GiNaC::matrix load_from_generate_cache(int integral, int block, time_t timestamp);
    void save_to_generate_cache(int integral, int block, time_t timestamp, const GiNaC::matrix& matrix);

    // subprocess management
    int max_subprocesses;
    int working_subprocesses;
    std::map<pid_t, std::pair<std::string, std::string>> read_subprocess_map;
    std::map<pid_t, std::string> expand_subprocess_map;
    std::map<pid_t, std::pair<int, int>> generate_subprocess_map;
    void read_subprocess_work(const std::string& key, const std::string& integral, const std::string& coefficient);
    void read_mainprocess_work(const std::string& key, const std::string& integral);
    void read_subprocess_yield(bool always_wait, const std::function<void(const std::string&, const std::string&)>& callback);
    void expand_subprocess_work(const std::string& key, const GiNaC::ex& coefficient, const GiNaC::lst& rules, int order, const GiNaC::symbol& eps);
    void expand_mainprocess_work(const std::string& key);
    void expand_subprocess_yield(bool always_wait);
    void generate_subprocess_work(int integral, int block, time_t timestamp, const GiNaC::matrix& raw_matrix, const GiNaC::lst& rules, const GiNaC::symbol& integral_symbol);
    void generate_mainprocess_work(int integral, int block, time_t timestamp, std::vector<std::vector<GiNaC::matrix>>* coefficient, std::vector<GiNaC::matrix>* bias);
    void generate_subprocess_yield(bool always_wait, time_t timestamp, std::vector<std::vector<GiNaC::matrix>>* coefficient, std::vector<GiNaC::matrix>* bias);

    GiNaC::ex get(GiNaC::symtab& _table, const std::string& _key,
                  const std::string& _prefix = "", 
                  const std::string& _suffix = "") {
        if (_table.find(_key) == _table.end()) {
            _table[_key] = GiNaC::symbol(_prefix + _key + _suffix);
        }
        return _table[_key];
    }

    std::string int_to_id(const std::string& _integral) {
        std::size_t start = _integral.find('[') + 1,
                    end = _integral.find(']');
        return _integral.substr(start, end - start);
    }

    void put_raw_symbols(const std::vector<std::string>& _keys) {
        for (auto& key: _keys)
            get(symbol_table, key);
    }

    GiNaC::symtab poly_to_terms(const std::string& base_indices, const GiNaC::ex& polynomial, GiNaC::symtab* record) {
        GiNaC::symtab result;
        auto indices = split(base_indices.c_str());
        auto termp = polynomial_iterator(polynomial), end = termp.end();
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
            result[repr] = factor;
            if (record != nullptr)
                (*record)[repr] = 0;
        }
        return result;
    }
};


#endif // CONFIG_HPP
