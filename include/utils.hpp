#ifndef UTILS_HPP
#define UTILS_HPP

#include <yaml-cpp/yaml.h>
#include <ginac/ginac.h>
#ifndef NO_GSL
#include <gsl/gsl_rng.h>
#endif // NO_GSL
#include <chrono>

#define START_TIME(description) auto description##_begin = std::chrono::high_resolution_clock::now()
#define END_TIME(description) auto description##_end = std::chrono::high_resolution_clock::now(); \
    std::chrono::duration<double, std::milli> description##_time_ms = description##_end - description##_begin
#define PRINT_TIME(description) std::cout << "Takes " << description##_time_ms.count() << " ms on " << #description << "()" << std::endl

std::vector<int> split(const char* str);
std::string combine(const std::vector<int>& values);

bool has_non_null_key(const YAML::Node& node, const std::string& key);

double to_double(const GiNaC::ex& ex);

GiNaC::matrix adjugate(const GiNaC::matrix& M);

#ifndef NO_GSL
// randomly generate Feynman parameters (x_1, ..., x_n)
// that are uniformly drawn from region 
// 0 <= x_i <= 1, x_1 + ... + x_n = 1
class random_feynman_params {
public:
    random_feynman_params(int n, unsigned long seed = gsl_rng_default_seed);
    ~random_feynman_params();
    void operator()(std::vector<double>& slots);
private:
    int length;
    std::vector<double> alpha;
    gsl_rng* rng;
};
#endif // NO_GSL

// correctly iterate over polynomial terms
class polynomial_iterator {
public:
    polynomial_iterator(const GiNaC::ex& polynomial);
    bool operator==(const polynomial_iterator& _other);
    bool operator!=(const polynomial_iterator& _other);
    GiNaC::ex operator*();
    polynomial_iterator& operator++();
    polynomial_iterator end();
private:
    const GiNaC::ex* polynomialp;
    bool has_multiple_terms;
    GiNaC::const_iterator iterator;
    GiNaC::const_iterator end_;
};

#endif // UTILS_HPP
