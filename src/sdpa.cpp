#include "sdpa.hpp"
#include "utils.hpp"

#ifndef NO_SDPA_LIB

sdpa_interface::sdpa_interface(const std::vector<std::vector<GiNaC::matrix>>& coefficients,
                               const std::vector<GiNaC::matrix>& bias,
                               const YAML::Node& config) {
    std::cerr << "SDPA start working..." << std::endl;
    fail = false;

    problem.setDisplay(nullptr);
    problem.setParameterType(SDPA::PARAMETER_DEFAULT);

    // set custom SDPA parameters
    if (has_non_null_key(config, "sdpa_params")) {
        auto sdpa_params = config["sdpa_params"];
        if (has_non_null_key(sdpa_params, "maxIteration")) {
            problem.setParameterMaxIteration(sdpa_params["maxIteration"].as<int>());
            std::cerr << "Set maxIteration = " << sdpa_params["maxIteration"].as<int>() << std::endl;
        }
        if (has_non_null_key(sdpa_params, "epsilonStar")) {
            problem.setParameterEpsilonStar(sdpa_params["epsilonStar"].as<double>());
            std::cerr << "Set epsilonStar = " << sdpa_params["epsilonStar"].as<double>() << std::endl;
        }
        if (has_non_null_key(sdpa_params, "epsilonDash")) {
            problem.setParameterEpsilonDash(sdpa_params["epsilonDash"].as<double>());
            std::cerr << "Set epsilonDash = " << sdpa_params["epsilonDash"].as<double>() << std::endl;
        }
        if (has_non_null_key(sdpa_params, "lambdaStar")) {
            problem.setParameterLambdaStar(sdpa_params["lambdaStar"].as<double>());
            std::cerr << "Set lambdaStar = " << sdpa_params["lambdaStar"].as<double>() << std::endl;
        }
        if (has_non_null_key(sdpa_params, "omegaStar")) {
            problem.setParameterOmegaStar(sdpa_params["omegaStar"].as<double>());
            std::cerr << "Set omegaStar = " << sdpa_params["omegaStar"].as<double>() << std::endl;
        }
        if (has_non_null_key(sdpa_params, "lowerBound")) {
            problem.setParameterLowerBound(sdpa_params["lowerBound"].as<double>());
            std::cerr << "Set lowerBound = " << sdpa_params["lowerBound"].as<double>() << std::endl;
        }
        if (has_non_null_key(sdpa_params, "upperBound")) {
            problem.setParameterUpperBound(sdpa_params["upperBound"].as<double>());
            std::cerr << "Set upperBound = " << sdpa_params["upperBound"].as<double>() << std::endl;
        }
        if (has_non_null_key(sdpa_params, "betaStar")) {
            problem.setParameterBetaStar(sdpa_params["betaStar"].as<double>());
            std::cerr << "Set betaStar = " << sdpa_params["betaStar"].as<double>() << std::endl;
        }
        if (has_non_null_key(sdpa_params, "betaBar")) {
            problem.setParameterBetaBar(sdpa_params["betaBar"].as<double>());
            std::cerr << "Set betaBar = " << sdpa_params["betaBar"].as<double>() << std::endl;
        }
        if (has_non_null_key(sdpa_params, "gammaStar")) {
            problem.setParameterGammaStar(sdpa_params["gammaStar"].as<double>());
            std::cerr << "Set gammaStar = " << sdpa_params["gammaStar"].as<double>() << std::endl;
        }
    }
    problem.printParameters(stdout);

    int nMasters = coefficients.size();
    int nBlock = bias.size();
    problem.inputConstraintNumber(nMasters + 1);
    problem.inputBlockNumber(nBlock);
    for (int i = 0; i < nBlock; i++) {
        problem.inputBlockSize(i + 1, bias[i].rows());
        problem.inputBlockType(i + 1, SDPA::SDP);
    }
    problem.initializeUpperTriangleSpace();

    for (int i = 0; i < nMasters; i++) {
        problem.inputCVec(i + 1, 0);
    }
    problem.inputCVec(nMasters + 1, 1);

    for (int j = 0; j < nBlock; j++) {
        int n = bias[j].rows();
        bias_matrices.push_back(std::vector<std::vector<GiNaC::ex>>(n, std::vector<GiNaC::ex>(n)));
        for (int k = 0; k < n; k++) {
            for (int l = k; l < n; l++) {
                auto element = bias[j](k, l);
                double value = 0.0;
                try {
                    value = -to_double(element);
                } catch (...) {
                    std::cerr << element << " is not a numeric value!" << std::endl;
                    fail = true;
                }
                problem.inputElement(0, j + 1, k + 1, l + 1, value);
                bias_matrices.back()[k][l] = element;
                bias_matrices.back()[l][k] = element;
            }
        }
    }

    for (int i = 0; i < nMasters; i++) {
        coefficient_matrices.push_back(std::vector<std::vector<std::vector<GiNaC::ex>>>(nBlock));
        for (int j = 0; j < nBlock; j++) {
            int n = coefficients[i][j].rows();
            coefficient_matrices.back()[j] = std::vector<std::vector<GiNaC::ex>>(n, std::vector<GiNaC::ex>(n));
            for (int k = 0; k < n; k++) {
                for (int l = k; l < n; l++) {
                    auto element = coefficients[i][j](k, l);
                    double value = 0.0;
                    try {
                        value = to_double(element);
                    } catch (...) {
                        std::cerr << element << " is not a numeric value!" << std::endl;
                        fail = true;
                    }
                    problem.inputElement(i + 1, j + 1, k + 1, l + 1, value);
                    coefficient_matrices.back()[j][k][l] = element;
                    coefficient_matrices.back()[j][l][k] = element;
                }
            }
        }
    }

    for (int j = 0; j < nBlock; j++) {
        int n = bias[j].rows();
        for (int k = 0; k < n; k++) {
            problem.inputElement(nMasters + 1, j + 1, k + 1, k + 1, 1);
        }
    }

    problem.initializeUpperTriangle();
    problem.initializeSolve();
}

void sdpa_interface::dump(std::ostream& out) {
    int nVariable = coefficient_matrices.size();
    int nBlock = bias_matrices.size();
    out << "b + Lambda";
    for (int i = 0; i < nVariable; i++)
        out << " + x" << i << " A" << i;
    out << " >= 0" << std::endl << std::endl;

    out << "b = {" << std::endl;
    for (int j = 0; j < nBlock; j++) {
        int n = bias_matrices[j].size();
        out << "    {{";
        for (int k = 0; k < n; k++) {
            for (int l = 0; l < n; l++) {
                out << bias_matrices[j][k][l];
                if (l != n - 1)
                    out << ", ";
            }
            if (k != n - 1)
                out << "},\n     {";
            else
                out << "}}," << std::endl;
        }
        out << std::endl;
    }
    out << "}" << std::endl << std::endl;

    for (int i = 0; i < nVariable; i++) {
        out << "A" << i << " = {" << std::endl;
        for (int j = 0; j < nBlock; j++) {
            int n = coefficient_matrices[i][j].size();
            out << "    {{";
            for (int k = 0; k < n; k++) {
                for (int l = 0; l < n; l++) {
                    out << coefficient_matrices[i][j][k][l];
                    if (l != n - 1)
                        out << ", ";
                }
                if (k != n - 1)
                    out << "},\n     {";
                else
                    out << "}}," << std::endl;
            }
            out << std::endl;
        }
        out << "}" << std::endl << std::endl;
    }
}

void sdpa_interface::solve() {
    if (fail) {
        std::cerr << "SDPA has failed during initialization!" << std::endl;
        std::cerr << "SDPA exiting..." << std::endl;
        return;
    }
    problem.solve();
    problem.printComputationTime(stdout);

    // write result
    result.stop_iteration = problem.getIteration();
    result.phase = problem.getPhaseValue();
    result.primal_obj = problem.getPrimalObj();
    result.dual_obj = problem.getDualObj();
    result.primal_err = problem.getPrimalError();
    result.dual_err = problem.getDualError();
    int n = problem.getConstraintNumber();
    result.x_vec = std::vector<double>(n);
    auto raw_x_vec = problem.getResultXVec();
    for (int i = 0; i < n; i++) {
        result.x_vec[i] = raw_x_vec[i];
    }
}

sdpa_interface::~sdpa_interface() {
    problem.terminate();
}

#else

#include <fstream>
#include <filesystem>

sdpa_interface::sdpa_interface(const std::vector<std::vector<GiNaC::matrix>>& coefficients,
                               const std::vector<GiNaC::matrix>& bias,
                               const YAML::Node& config) {
    std::cerr << "SDPA start working..." << std::endl;
    fail = false;


    std::filesystem::create_directory("logs");
    std::ofstream param_file(std::filesystem::path("logs").append("param.sdpa"));

    // set custom SDPA parameters
    if (has_non_null_key(config, "sdpa_params")) {
        auto sdpa_params = config["sdpa_params"];
        if (has_non_null_key(sdpa_params, "maxIteration")) {
            param_file << sdpa_params["maxIteration"].as<int>() << "\tunsigned int maxIteration;" << std::endl;
            std::cerr << "Set maxIteration = " << sdpa_params["maxIteration"].as<int>() << std::endl;
        } else {
            param_file << "100\tunsigned int maxIteration;" << std::endl;
        }
        if (has_non_null_key(sdpa_params, "epsilonStar")) {
            param_file << sdpa_params["epsilonStar"].as<double>() << "\tdouble 0.0 < epsilonStar;" << std::endl;
            std::cerr << "Set epsilonStar = " << sdpa_params["epsilonStar"].as<double>() << std::endl;
        } else {
            param_file << "1.0E-7\tdouble 0.0 < epsilonStar;" << std::endl;
        }
        if (has_non_null_key(sdpa_params, "lambdaStar")) {
            param_file << sdpa_params["lambdaStar"].as<double>() << "\tdouble 0.0 < lambdaStar;" << std::endl;
            std::cerr << "Set lambdaStar = " << sdpa_params["lambdaStar"].as<double>() << std::endl;
        } else {
            param_file << "1.0E2\tdouble 0.0 < lambdaStar;" << std::endl;
        }
        if (has_non_null_key(sdpa_params, "omegaStar")) {
            param_file << sdpa_params["omegaStar"].as<double>() << "\tdouble 1.0 < omegaStar;" << std::endl;
            std::cerr << "Set omegaStar = " << sdpa_params["omegaStar"].as<double>() << std::endl;
        } else {
            param_file << "2.0\tdouble 1.0 < omegaStar;" << std::endl;
        }
        if (has_non_null_key(sdpa_params, "lowerBound")) {
            param_file << sdpa_params["lowerBound"].as<double>() << "\tdouble lowerBound;" << std::endl;
            std::cerr << "Set lowerBound = " << sdpa_params["lowerBound"].as<double>() << std::endl;
        } else {
            param_file << "-1.0E5\tdouble lowerBound;" << std::endl;
        }
        if (has_non_null_key(sdpa_params, "upperBound")) {
            param_file << sdpa_params["upperBound"].as<double>() << "\tdouble upperBound;" << std::endl;
            std::cerr << "Set upperBound = " << sdpa_params["upperBound"].as<double>() << std::endl;
        } else {
            param_file << "1.0E5\tdouble upperBound;" << std::endl;
        }
        if (has_non_null_key(sdpa_params, "betaStar")) {
            param_file << sdpa_params["betaStar"].as<double>() << "\tdouble 0.0 <= betaStar <  1.0;" << std::endl;
            std::cerr << "Set betaStar = " << sdpa_params["betaStar"].as<double>() << std::endl;
        } else {
            param_file << "0.1\tdouble 0.0 <= betaStar <  1.0;" << std::endl;
        }
        if (has_non_null_key(sdpa_params, "betaBar")) {
            param_file << sdpa_params["betaBar"].as<double>() << "\tdouble 0.0 <= betaBar  <  1.0, betaStar <= betaBar;" << std::endl;
            std::cerr << "Set betaBar = " << sdpa_params["betaBar"].as<double>() << std::endl;
        } else {
            param_file << "0.2\tdouble 0.0 <= betaBar  <  1.0, betaStar <= betaBar;" << std::endl;
        }
        if (has_non_null_key(sdpa_params, "gammaStar")) {
            param_file << sdpa_params["gammaStar"].as<double>() << "\tdouble 0.0 < gammaStar  <  1.0;" << std::endl;
            std::cerr << "Set gammaStar = " << sdpa_params["gammaStar"].as<double>() << std::endl;
        } else {
            param_file << "0.9\tdouble 0.0 < gammaStar  <  1.0;" << std::endl;
        }
        if (has_non_null_key(sdpa_params, "epsilonDash")) {
            param_file << sdpa_params["epsilonDash"].as<double>() << "\tdouble 0.0 < epsilonDash;" << std::endl;
            std::cerr << "Set epsilonDash = " << sdpa_params["epsilonDash"].as<double>() << std::endl;
        } else {
            param_file << "1.0E-7\tdouble 0.0 < epsilonDash;" << std::endl;
        }
        if (has_non_null_key(sdpa_params, "xPrint")) {
            param_file << sdpa_params["xPrint"].as<std::string>() << "\tchar* xPrint\t(default %+8.3e,   NOPRINT skips printout)" << std::endl;
            std::cerr << "Set xPrint = " << sdpa_params["xPrint"].as<std::string>() << std::endl;
        } else {
            param_file << "%+8.3e\tchar* xPrint\t(default %+8.3e,   NOPRINT skips printout)" << std::endl;
        }
        if (has_non_null_key(sdpa_params, "XPrint")) {
            param_file << sdpa_params["XPrint"].as<std::string>() << "\tchar* XPrint\t(default %+8.3e,   NOPRINT skips printout)" << std::endl;
            std::cerr << "Set XPrint = " << sdpa_params["XPrint"].as<std::string>() << std::endl;
        } else {
            param_file << "%+8.3e\tchar* XPrint\t(default %+8.3e,   NOPRINT skips printout)" << std::endl;
        }
        if (has_non_null_key(sdpa_params, "YPrint")) {
            param_file << sdpa_params["YPrint"].as<std::string>() << "\tchar* YPrint\t(default %+8.3e,   NOPRINT skips printout)" << std::endl;
            std::cerr << "Set YPrint = " << sdpa_params["YPrint"].as<std::string>() << std::endl;
        } else {
            param_file << "%+8.3e\tchar* YPrint\t(default %+8.3e,   NOPRINT skips printout)" << std::endl;
        }
        if (has_non_null_key(sdpa_params, "infPrint")) {
            param_file << sdpa_params["infPrint"].as<std::string>() << "\tchar* infPrint\t(default %+10.16e, NOPRINT skips printout)" << std::endl;
            std::cerr << "Set infPrint = " << sdpa_params["infPrint"].as<std::string>() << std::endl;
        } else {
            param_file << "%+10.16e\tchar* infPrint\t(default %+10.16e, NOPRINT skips printout)" << std::endl;
        }
    }
    param_file.close();

    std::ofstream problem_file(std::filesystem::path("logs").append("problem.in"));
    int nMasters = coefficients.size();
    int nBlock = bias.size();
    problem_file << "    " << (nMasters + 1) << " = mDIM" << std::endl;
    problem_file << "    " << nBlock << " = nBLOCK" << std::endl;
    problem_file << "    ";
    for (int i = 0; i < nBlock; i++)
        problem_file << bias[i].rows() << "    ";
    problem_file << " = bLOCKsTRUCT" << std::endl;

    problem_file << "{";
    for (int i = 0; i < nMasters; i++) {
        problem_file << "0, ";
    }
    problem_file << "1}" << std::endl;

    problem_file << "{" << std::endl;
    for (int j = 0; j < nBlock; j++) {
        int n = bias[j].rows();
        bias_matrices.push_back(std::vector<std::vector<GiNaC::ex>>(n, std::vector<GiNaC::ex>(n)));
        problem_file << "{ ";
        for (int k = 0; k < n; k++) {
            if (k == 0)
                problem_file << "{";
            else 
                problem_file << "},\n  {";
            for (int l = 0; l < n; l++) {
                auto element = bias[j](k, l);                
                double value = 0.0;
                try {
                    value = -to_double(element);
                    problem_file << GiNaC::ex_to<GiNaC::numeric>(-element).evalf();
                } catch (...) {
                    std::cerr << element << " is not a numeric value!" << std::endl;
                    fail = true;
                }
                if (l != n - 1)
                    problem_file << ", ";
                bias_matrices.back()[k][l] = element;
                bias_matrices.back()[l][k] = element;
            }
        }
        problem_file << "} }" << std::endl;
    }
    problem_file << "}" << std::endl;

    for (int i = 0; i < nMasters; i++) {
        coefficient_matrices.push_back(std::vector<std::vector<std::vector<GiNaC::ex>>>(nBlock));
        problem_file << "{" << std::endl;
        for (int j = 0; j < nBlock; j++) {
            int n = coefficients[i][j].rows();
            coefficient_matrices.back()[j] = std::vector<std::vector<GiNaC::ex>>(n, std::vector<GiNaC::ex>(n));
            problem_file << "{ ";
            for (int k = 0; k < n; k++) {
                if (k == 0)
                    problem_file << "{";
                else 
                    problem_file << "},\n  {";
                for (int l = 0; l < n; l++) {
                    auto element = coefficients[i][j](k, l);
                    double value = 0.0;
                    try {
                        value = to_double(element);
                        problem_file << GiNaC::ex_to<GiNaC::numeric>(element).evalf();
                    } catch (...) {
                        std::cerr << element << " is not a numeric value!" << std::endl;
                        fail = true;
                    }
                    if (l != n - 1)
                        problem_file << ", ";
                    coefficient_matrices.back()[j][k][l] = element;
                    coefficient_matrices.back()[j][l][k] = element;
                }
            }
            problem_file << "} }" << std::endl;
        }
        problem_file << "}" << std::endl;
    }

    problem_file << "{" << std::endl;
    for (int j = 0; j < nBlock; j++) {
        int n = bias[j].rows();
        problem_file << "{ ";
        for (int k = 0; k < n; k++) {
            if (k == 0)
                problem_file << "{";
            else 
                problem_file << "},\n  {";
            for (int l = 0; l < n; l++) {
                if (l == k)
                    problem_file << 1;
                else
                    problem_file << 0;
                if (l != n - 1)
                    problem_file << ", ";
            }
        }
        problem_file << "} }" << std::endl;
    }
    problem_file << "}" << std::endl;
    
    problem_file.close();
}

void sdpa_interface::dump(std::ostream& out) {
    int nVariable = coefficient_matrices.size();
    int nBlock = bias_matrices.size();
    out << "b + Lambda";
    for (int i = 0; i < nVariable; i++)
        out << " + x" << i << " A" << i;
    out << " >= 0" << std::endl << std::endl;

    out << "b = {" << std::endl;
    for (int j = 0; j < nBlock; j++) {
        int n = bias_matrices[j].size();
        out << "    {{";
        for (int k = 0; k < n; k++) {
            for (int l = 0; l < n; l++) {
                out << bias_matrices[j][k][l];
                if (l != n - 1)
                    out << ", ";
            }
            if (k != n - 1)
                out << "},\n     {";
            else
                out << "}}," << std::endl;
        }
        out << std::endl;
    }
    out << "}" << std::endl << std::endl;

    for (int i = 0; i < nVariable; i++) {
        out << "A" << i << " = {" << std::endl;
        for (int j = 0; j < nBlock; j++) {
            int n = coefficient_matrices[i][j].size();
            out << "    {{";
            for (int k = 0; k < n; k++) {
                for (int l = 0; l < n; l++) {
                    out << coefficient_matrices[i][j][k][l];
                    if (l != n - 1)
                        out << ", ";
                }
                if (k != n - 1)
                    out << "},\n     {";
                else
                    out << "}}," << std::endl;
            }
            out << std::endl;
        }
        out << "}" << std::endl << std::endl;
    }
}

void sdpa_interface::solve() {
}

sdpa_interface::~sdpa_interface() {
}

#endif // NO_SDPA_LIB
