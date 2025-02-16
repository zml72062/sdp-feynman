#include "sdpa.hpp"
#include "utils.hpp"

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
            }
        }
    }

    for (int i = 0; i < nMasters; i++) {
        for (int j = 0; j < nBlock; j++) {
            int n = coefficients[i][j].rows();
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
