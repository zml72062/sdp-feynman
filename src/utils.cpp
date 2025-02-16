#include <gsl/gsl_randist.h>
#include "utils.hpp"

std::vector<int> split(const char* str) {
    std::vector<int> nums;
    const char* ptr = str;
    char* endptr;
    while (true) {
        nums.push_back(std::strtol(ptr, &endptr, 10));
        if (*endptr == '\0')
            return nums;
        ptr = endptr + 1;
    }
}

std::string combine(const std::vector<int>& values) {
    std::string str = "";
    int num_values = values.size();
    for (int i = 0; i < num_values; i++) {
        str += std::to_string(values[i]);
        if (i != num_values - 1)
            str += ",";
    }
    return str;
}

bool has_non_null_key(const YAML::Node& node, const std::string& key) {
    auto node_map = node.as<std::map<std::string, YAML::Node>>();
    return node_map.find(key) != node_map.end() && node[key].Type() != YAML::NodeType::Null;
}

double to_double(const GiNaC::ex& ex) {
    return GiNaC::ex_to<GiNaC::numeric>(GiNaC::ex_to<GiNaC::numeric>(ex).evalf()).to_double();
}

GiNaC::matrix adjugate(const GiNaC::matrix& M) {
    auto n = M.rows(), m = M.cols();
    if (n != m)
        throw std::runtime_error("adjugate(): non-square matrix");
    
    GiNaC::matrix adj(n, n);
    if (n == 1) {
        adj(0, 0) = 1;
        return adj;
    }

    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = 0; j < n; j++) {
            adj(i, j) = GiNaC::pow(-1, i + j) * 
                GiNaC::ex_to<GiNaC::matrix>(
                    GiNaC::reduced_matrix(M, j, i)
                ).determinant();
        }
    }

    return adj;
}

random_feynman_params::random_feynman_params(int n, unsigned long seed) : length(n), alpha(n, 1.0) {
    gsl_rng_env_setup();
    rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, seed);
}

random_feynman_params::~random_feynman_params() {
    gsl_rng_free(rng);
}

void random_feynman_params::operator()(std::vector<double>& slots) {
    gsl_ran_dirichlet(rng, length, alpha.data(), slots.data());
}

polynomial_iterator::polynomial_iterator(const GiNaC::ex& polynomial)
    : polynomialp(&polynomial) {
    
    has_multiple_terms = GiNaC::is_exactly_a<GiNaC::add>(polynomial);
    if (has_multiple_terms) // need iteration
        iterator = polynomial.begin();
    else // do not need iteration
        iterator = polynomial.end() - 1;
    end_ = polynomial.end();
}

bool polynomial_iterator::operator==(const polynomial_iterator& _other) {
    return iterator == _other.iterator;
}

bool polynomial_iterator::operator!=(const polynomial_iterator& _other) {
    return iterator != _other.iterator;
}

GiNaC::ex polynomial_iterator::operator*() {
    if (has_multiple_terms)
        return *iterator;
    return *polynomialp;
}

polynomial_iterator& polynomial_iterator::operator++() {
    ++iterator;
    return *this;
}

polynomial_iterator polynomial_iterator::end() {
    auto end_iterator = polynomial_iterator(*polynomialp);
    end_iterator.iterator = polynomialp->end();
    return end_iterator;
}
