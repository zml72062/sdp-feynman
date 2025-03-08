#ifndef SDPA_HPP
#define SDPA_HPP

#include <iostream>
#ifndef NO_SDPA_LIB
#include <cstdio>
#include <cstdlib>
#include <sdpa_call.h>
#endif // NO_SDPA_LIB
#include <ginac/ginac.h>
#include <yaml-cpp/yaml.h>

#ifndef NO_SDPA_LIB
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
    void dump(std::ostream& out);
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
    // for debug
    std::vector<std::vector<std::vector<std::vector<GiNaC::ex>>>> coefficient_matrices;
    std::vector<std::vector<std::vector<GiNaC::ex>>> bias_matrices;
};

#else

class sdpa_interface {
public:
    sdpa_interface(const std::vector<std::vector<GiNaC::matrix>>& coefficients,
                   const std::vector<GiNaC::matrix>& bias,
                   const YAML::Node& config);
    void solve();
    void dump(std::ostream& out);
    ~sdpa_interface();

    bool get_fail() {
        return fail;
    }
private:
    // flag indicating whether error has occurred during initialization
    bool fail;
    // for debug
    std::vector<std::vector<std::vector<std::vector<GiNaC::ex>>>> coefficient_matrices;
    std::vector<std::vector<std::vector<GiNaC::ex>>> bias_matrices;
};

#endif // NO_SDPA_LIB

#endif // SDPA_HPP
