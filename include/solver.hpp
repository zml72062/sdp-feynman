#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <ginac/ginac.h>
#include <yaml-cpp/yaml.h>

class master_solver {
public:
    master_solver(const std::vector<std::string>& effective_master_table,
                  std::vector<GiNaC::symtab>& numeric_integral_table,
                  const YAML::Node& config, bool dump)
        : effective_master_tablep(&effective_master_table),
          numeric_integral_tablep(&numeric_integral_table),
          configp(&config), will_dump(dump) { 
        
        fail = false;
        int order = numeric_integral_table.size();
        for (int i = 0; i < order; i++) {
            for (auto& name: effective_master_table) {
                variables_to_solve.append(numeric_integral_table[i][name]);
            }
        }
    }
    
    void solve_from(const std::vector<GiNaC::matrix>& matrices, void* config_parserp);
    
#ifndef NO_SDPA_LIB
    const GiNaC::lst& get_result() {
        return computed_values;
    }

    bool get_fail() {
        return fail;
    }
#endif // NO_SDPA_LIB
private:
    const std::vector<std::string>* effective_master_tablep;
    std::vector<GiNaC::symtab>* numeric_integral_tablep;
    const YAML::Node* configp;

    GiNaC::lst variables_to_solve;
    GiNaC::lst computed_values;
    bool fail;
    bool will_dump;
    const double positivity_threshold = 1e-5;
};

#endif // SOLVER_HPP

