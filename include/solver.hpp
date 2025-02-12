#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <ginac/ginac.h>
#include <yaml-cpp/yaml.h>

class master_solver {
public:
    master_solver(const std::vector<std::string>& effective_master_table,
                  std::vector<GiNaC::symtab>& numeric_integral_table,
                  const YAML::Node& config)
        : effective_master_tablep(&effective_master_table),
          numeric_integral_tablep(&numeric_integral_table),
          configp(&config) { 
        
        fail = false;
        int order = numeric_integral_table.size();
        for (int i = 0; i < order; i++) {
            for (auto& name: effective_master_table) {
                variables_to_solve.append(numeric_integral_table[i][name]);
            }
        }
    }
    
    void solve_from(const std::vector<GiNaC::matrix>& matrices);
    
    const GiNaC::lst& get_result() {
        return computed_values;
    }

    bool get_fail() {
        return fail;
    }
private:
    const std::vector<std::string>* effective_master_tablep;
    std::vector<GiNaC::symtab>* numeric_integral_tablep;
    const YAML::Node* configp;

    GiNaC::lst variables_to_solve;
    GiNaC::lst computed_values;
    bool fail;
    const double positivity_threshold = 1e-5;
};

#endif

