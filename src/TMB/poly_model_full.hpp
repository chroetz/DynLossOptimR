// poly_model_full.hpp
// Do NOT #include <TMB.hpp>

#include <algorithm>

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template <class Type>
void fill_poly_cols_recursive_full(
    matrix<Type>& featureMatrix,
    int& col_idx,
    const vector<Type>& current_col,
    int var_idx,
    int d,
    int remaining_deg,
    const matrix<Type>& inputMatrix
) {
    if (col_idx >= featureMatrix.cols()) return;
    featureMatrix.col(col_idx) = current_col;
    col_idx++;
    if (remaining_deg <= 0) return;
    for (int i = var_idx; i < d; i++) {
        vector<Type> next_col = (current_col.array() * inputMatrix.col(i).array()).matrix();
        fill_poly_cols_recursive_full(featureMatrix, col_idx, next_col, i, d, remaining_deg - 1, inputMatrix);
    }
}

template <class Type>
void fill_poly_full(matrix<Type>& featureMatrix, const matrix<Type>& inputMatrix, int deg) {
    int col_pos = 0;
    int n = inputMatrix.rows();
    int d = inputMatrix.cols();
    vector<Type> start_col(n);
    start_col.setOnes();
    fill_poly_cols_recursive_full(featureMatrix, col_pos, start_col, 0, d, deg, inputMatrix);
}

template<class Type>
Type poly_model_full(objective_function<Type>* obj) {
    DATA_MATRIX(obs);
    DATA_VECTOR(weights_obs);
    DATA_VECTOR(weights_ana);
    DATA_VECTOR(weights_pen);
    DATA_VECTOR(weights_coef);
    DATA_INTEGER(deg);

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
        fill_poly_full(features, currentInput, deg);
        currentInput = features * coef;
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
