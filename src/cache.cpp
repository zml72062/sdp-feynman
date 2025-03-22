#include "config.hpp"
#include "utils.hpp"
#include <fstream>
#include <filesystem>

bool config_parser::read_cache_exists(const std::string& key, const std::string& integral) {
    return std::filesystem::exists(std::filesystem::path(cache_dir).append("read").append("cache_" + key + "_" + integral));
}

GiNaC::ex config_parser::load_from_read_cache(const std::string& key, const std::string& integral) {
    GiNaC::lst syms;
    for (auto& symbol: symbol_table) {
        syms.append(symbol.second);
    }
    
    GiNaC::archive ar;
    std::ifstream in(std::filesystem::path(cache_dir).append("read").append("cache_" + key + "_" + integral),
                     std::ios::binary);
    in >> ar;
    in.close();
    return ar.unarchive_ex(syms, "coeff");
}

void config_parser::save_to_read_cache(const std::string& key, const std::string& integral, const GiNaC::ex& coefficient) {
    GiNaC::archive ar;
    ar.archive_ex(coefficient, "coeff");
    std::ofstream out(std::filesystem::path(cache_dir).append("read").append("cache_" + key + "_" + integral),
                      std::ios::binary);
    out << ar;
    out.close();
}

GiNaC::ex config_parser::read_ibp_simple(const std::string& key, const std::string& integral) {
    auto coefficient = load_from_read_cache(key, integral);
    auto key_indices = split(key.c_str());
    auto integral_indices = split(integral.c_str());
    int n_indices = key_indices.size();
    int sum_diff = 0;
    for (int i = 0; i < n_indices; i++) {
        sum_diff += (key_indices[i] - integral_indices[i]);
    }
    return coefficient * GiNaC::pow(-1, sum_diff);
}

bool config_parser::expand_cache_exists(const std::string& key) {
    return std::filesystem::exists(std::filesystem::path(cache_dir).append("expand").append("cache_" + key));
}


GiNaC::ex config_parser::load_from_expand_cache(const std::string& key) {
    GiNaC::lst syms;
    int order = numeric_integral_table.size();
    for (int i = 0; i < order; i++) {
        for (auto& name: effective_master_table) {
            syms.append(numeric_integral_table[i][name]);
        }
    }
    
    GiNaC::archive ar;
    std::ifstream in(std::filesystem::path(cache_dir).append("expand").append("cache_" + key),
                     std::ios::binary);
    in >> ar;
    in.close();
    return ar.unarchive_ex(syms, "coeff");
}

void config_parser::save_to_expand_cache(const std::string& key, const GiNaC::ex& coefficient) {
    GiNaC::archive ar;
    ar.archive_ex(coefficient, "coeff");
    std::ofstream out(std::filesystem::path(cache_dir).append("expand").append("cache_" + key),
                      std::ios::binary);
    out << ar;
    out.close();
}

GiNaC::matrix config_parser::load_from_generate_cache(int integral, int block, time_t timestamp) {
    GiNaC::lst syms;
    GiNaC::archive ar;
    std::ifstream in(std::filesystem::path(cache_dir).append("generate")
                     .append("cache_" + std::to_string(integral) + "_" + std::to_string(block) + "_" + std::to_string(timestamp)),
                     std::ios::binary);
    in >> ar;
    in.close();
    return GiNaC::ex_to<GiNaC::matrix>(ar.unarchive_ex(syms, "coeff"));
}

void config_parser::save_to_generate_cache(int integral, int block, time_t timestamp, const GiNaC::matrix& matrix) {
    GiNaC::archive ar;
    ar.archive_ex(matrix, "coeff");
    std::ofstream out(std::filesystem::path(cache_dir).append("generate")
                      .append("cache_" + std::to_string(integral) + "_" + std::to_string(block) + "_" + std::to_string(timestamp)),
                      std::ios::binary);
    out << ar;
    out.close();
}
