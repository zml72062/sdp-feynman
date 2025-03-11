#include "config.hpp"
#include <fstream>
#include <filesystem>

bool config_parser::cache_exists(const std::string& key, const std::string& integral) {
    return std::filesystem::exists(std::filesystem::path(cache_dir).append("cache_" + key + "_" + integral));
}

GiNaC::ex config_parser::load_from_cache(const std::string& key, const std::string& integral) {
    GiNaC::lst syms;
    for (auto& symbol: symbol_table) {
        syms.append(symbol.second);
    }
    
    GiNaC::archive ar;
    std::ifstream in(std::filesystem::path(cache_dir).append("cache_" + key + "_" + integral),
                     std::ios::binary);
    in >> ar;
    in.close();
    return ar.unarchive_ex(syms, "coeff");
}

void config_parser::save_to_cache(const std::string& key, const std::string& integral, const GiNaC::ex& coefficient) {
    GiNaC::archive ar;
    ar.archive_ex(coefficient, "coeff");
    std::ofstream out(std::filesystem::path(cache_dir).append("cache_" + key + "_" + integral),
                      std::ios::binary);
    out << ar;
    out.close();
}
