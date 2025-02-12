#ifndef SDPA_HPP
#define SDPA_HPP

#include <cstdio>
#include <cstdlib>
#include <sdpa_call.h>
#include <ginac/ginac.h>
#include <yaml-cpp/yaml.h>

struct sdpa_result {
    int stop_iteration;
    SDPA::PhaseType phase;
    double primal_obj;
    double dual_obj;
    double primal_err;
    double dual_err;
    std::vector<double> x_vec;
};

class sdpa_interface {
public:
    sdpa_interface(const std::vector<std::vector<GiNaC::matrix>>& coefficients,
                   const std::vector<GiNaC::matrix>& bias,
                   const YAML::Node& config);
    void solve();
    ~sdpa_interface();

    const sdpa_result& get_result() {
        return result;
    }

    bool get_fail() {
        return fail;
    }
private:
    // data structure for the SDPA problem
    SDPA problem;
    // data structure for the solution
    sdpa_result result;
    // flag indicating whether error has occurred during initialization
    bool fail;
};


#endif
