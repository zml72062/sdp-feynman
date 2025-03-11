#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <yaml-cpp/yaml.h>
#include <ginac/ginac.h>
#include <sys/types.h>
#include <unistd.h>
#include "parse.hpp"
#include "solver.hpp"


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
    bool cache_exists(const std::string& key, const std::string& integral);
    GiNaC::ex load_from_cache(const std::string& key, const std::string& integral);
    void save_to_cache(const std::string& key, const std::string& integral, const GiNaC::ex& coefficient);

    // subprocess management
    int max_subprocesses;
    int working_subprocesses;
    std::map<pid_t, std::pair<std::string, std::string>> subprocess_map;
    void subprocess_work(const std::string& key, const std::string& integral, const std::string& coefficient);
    void mainprocess_work(const std::string& key, const std::string& integral);
    void subprocess_yield(bool always_wait);

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
};


#endif // CONFIG_HPP
