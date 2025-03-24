#include "asy.hpp"
#include "utils.hpp"
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <fstream>
#include <filesystem>


GiNaC::matrix asy::get_term_orders(const GiNaC::ex& polynomial) {
    std::vector<std::vector<int>> result;
    auto termp = polynomial_iterator(polynomial), end = termp.end();
    for (; termp != end; ++termp) {
        result.push_back(std::vector<int>());
        auto term = *termp;
        for (auto& i: *effective_feynman_paramsp) {
            result.back().push_back(term.degree((*feynman_paramsp)[i]));
            term = term.lcoeff((*feynman_paramsp)[i]);
        }
        result.back().push_back(term.degree(*xp));
    }

    if (result.size() == 0)
        return GiNaC::matrix(0, 0);
    
    int row = result.size(), col = result[0].size();
    GiNaC::matrix matrix(row, col);
    for (int i = 0; i < row; i++)
        for (int j = 0; j < col; j++)
            matrix(i, j) = result[i][j];
    return matrix;
}


GiNaC::matrix asy::export_to_python() {
    return get_term_orders((U_poly * F_poly).expand());
}


void asy::python_work(char* python_path) {
    pid_t pid = fork();

    START_TIME(asy);
    if (pid == 0) { // child
        std::filesystem::create_directory("logs");
        int fd = open(std::filesystem::path("logs").append("asy.input").c_str(), 
                      O_RDWR | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);
        dup2(fd, STDOUT_FILENO);
        close(fd);
        std::cout << export_to_python() << std::endl;
        fd = open(std::filesystem::path("logs").append("asy.output").c_str(), 
                  O_RDWR | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);
        dup2(fd, STDOUT_FILENO);
        close(fd);
        char* const argv[4] = {python_path, (char*)"asy.py", 
                               (char*)std::filesystem::path("logs").append("asy.input").c_str(), 0};
        if (execv(python_path, argv) < 0) {
            std::cerr << "Error occurred while executing Python!" << std::endl;
        }
        exit(0);
    } else { // parent
        waitpid(pid, 0, 0);
        std::ifstream asy_result_file(std::filesystem::path("logs").append("asy.output"));
        int row = 0, col = 0;
        asy_result_file >> row >> col;
        if (row > 0 && col > 0) {
            scaling_vectors = GiNaC::matrix(row, col);
            GiNaC::parser parser;
            std::string input;
            for (int i = 0; i < row; i++) {
                for (int j = 0; j < col; j++) {
                    asy_result_file >> input;
                    scaling_vectors(i, j) = parser(input);
                }
            }
        }
        asy_result_file.close();
        END_TIME(asy);
        PRINT_TIME(asy);
    }
}


void asy::compute_asymptotic_polys() {
    auto U_term_orders = get_term_orders(U_poly);
    auto F_term_orders = get_term_orders(F_poly);
    int U_terms = U_term_orders.rows();
    int F_terms = F_term_orders.rows();

    int num_vectors = scaling_vectors.rows(), n = scaling_vectors.cols();
    std::vector<double> U_term_power(U_terms), F_term_power(F_terms);
    double U_term_min_power, F_term_min_power;
    for (int i = 0; i < num_vectors; i++) {
        U_term_min_power = 1e10;
        // find minimum power in U polynomial
        for (int j = 0; j < U_terms; j++) {
            double sum = 0;
            for (int k = 0; k < n; k++) {
                sum += to_double(U_term_orders(j, k) * scaling_vectors(i, k));
            }
            sum += to_double(U_term_orders(j, n));
            U_term_power[j] = sum;
            if (sum < U_term_min_power)
                U_term_min_power = sum;
        }

        auto termp = polynomial_iterator(U_poly), end = termp.end();
        int j = 0;
        GiNaC::ex asymptotic_poly = 0;
        for (; termp != end; ++termp, ++j) {
            auto term = *termp;
            if (U_term_power[j] > U_term_min_power - 1e-8 &&
                U_term_power[j] < U_term_min_power + 1e-8)
                asymptotic_poly += term;
        }
        asymptotic_Us.append(asymptotic_poly);
    }

    for (int i = 0; i < num_vectors; i++) {
        F_term_min_power = 1e10;
        // find minimum power in F polynomial
        for (int j = 0; j < F_terms; j++) {
            double sum = 0;
            for (int k = 0; k < n; k++) {
                sum += to_double(F_term_orders(j, k) * scaling_vectors(i, k));
            }
            sum += to_double(F_term_orders(j, n));
            F_term_power[j] = sum;
            if (sum < F_term_min_power)
                F_term_min_power = sum;
        }

        auto termp = polynomial_iterator(F_poly), end = termp.end();
        int j = 0;
        GiNaC::ex asymptotic_poly = 0;
        for (; termp != end; ++termp, ++j) {
            auto term = *termp;
            if (F_term_power[j] > F_term_min_power - 1e-8 &&
                F_term_power[j] < F_term_min_power + 1e-8)
                asymptotic_poly += term;
        }
        asymptotic_Fs.append(asymptotic_poly);
    }
}

