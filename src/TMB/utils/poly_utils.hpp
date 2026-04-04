// poly_utils.hpp
// Shared polynomial feature utilities used by poly_model.hpp and poly_model_full.hpp.
// Do NOT #include <TMB.hpp>

#ifndef POLY_UTILS_HPP
#define POLY_UTILS_HPP

template <class Type>
void fill_poly_cols_recursive(
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
        fill_poly_cols_recursive(featureMatrix, col_idx, next_col, i, d, remaining_deg - 1, inputMatrix);
    }
}

template <class Type>
void fill_poly(matrix<Type>& featureMatrix, const matrix<Type>& inputMatrix, int deg) {
    int col_pos = 0;
    int n = inputMatrix.rows();
    int d = inputMatrix.cols();
    vector<Type> start_col(n);
    start_col.setOnes();
    fill_poly_cols_recursive(featureMatrix, col_pos, start_col, 0, d, deg, inputMatrix);
}

#endif // POLY_UTILS_HPP
