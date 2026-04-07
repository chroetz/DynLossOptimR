// poly_model_full.hpp
// Do NOT #include <TMB.hpp>

#include <algorithm>
#include "utils/poly_utils.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type poly_model_full(objective_function<Type>* obj) {
    DATA_MATRIX(obs);
    DATA_VECTOR(weights_obs);
    DATA_VECTOR(weights_ana);
    DATA_VECTOR(weights_pen);
    DATA_VECTOR(weights_coef);
    DATA_INTEGER(deg);
    DATA_INTEGER(intermediate);

    PARAMETER_MATRIX(ana);
    PARAMETER_MATRIX(coef);

    int n = obs.rows();
    int p = coef.rows();
    int layers = std::max((int)weights_obs.size() - 1, (int)weights_ana.size());
    int penalty_order = weights_pen.size();

    matrix<Type> diff0 = ana - obs;
    Type nll = weights_obs[0] * (diff0.array() * diff0.array()).sum() / n;

    matrix<Type> currentInput = ana;
    matrix<Type> features(n, p);

    for(int k = 0; k < layers; k++) {
        for (int j = 0; j < intermediate; j++) {
            fill_poly(features, currentInput, deg);
            currentInput = features * coef;
        }
        matrix<Type> diff = currentInput.topRows(n-k-1) - obs.bottomRows(n-k-1);
        if (weights_obs.size() > k+1) {
            nll += weights_obs[k+1] * (diff.array() * diff.array()).sum() / (n-k-1);
        }
        diff = currentInput.topRows(n-k-1) - ana.bottomRows(n-k-1);
        if (weights_ana.size() > k) {
            nll += weights_ana[k] * (diff.array() * diff.array()).sum() / (n-k-1);
        }
    }

    matrix<Type> current = ana;
    for(int j = 0; j < penalty_order; j++) {
        int m = current.rows();
        current = current.topRows(m-1) - current.bottomRows(m-1);
        nll += weights_pen[j] * (current.array() * current.array()).sum() / (m-1);
    }

    for(int j = 0; j < p; j++) {
        nll += weights_coef[j] * coef.row(j).dot(coef.row(j)) / p;
    }

    return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
