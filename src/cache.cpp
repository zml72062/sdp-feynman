#include "config.hpp"
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
