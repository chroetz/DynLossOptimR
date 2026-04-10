// poly_model_coef.hpp
// Do NOT #include <TMB.hpp>

#include "utils/poly_utils.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type poly_model_coef0(objective_function<Type>* obj) {
    DATA_MATRIX(obs);
    DATA_INTEGER(deg);
    DATA_INTEGER(layers);

    PARAMETER_MATRIX(coef);

    Type scale = Type(1e4); // to soft squash state to [-1e4, 1e4]

    int n = obs.rows();
    int p = coef.rows();

    Type nll = 0;

    matrix<Type> currentInput = obs;
    matrix<Type> features(n, p);

    for(int k = 0; k < layers; k++) {
        fill_poly(features, currentInput, deg);
        currentInput = features * coef;
        currentInput = scale * (currentInput / scale).array().tanh(); // Soft squash to [-scale, scale].
        matrix<Type> diff = currentInput.topRows(n-k-1) - obs.bottomRows(n-k-1);
        nll += (diff.array() * diff.array()).sum() / (n-k-1);
    }

    return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
