// src/TMB/poly_model.hpp
// Do NOT #include <TMB.hpp>

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template <class Type>
void fill_poly_cols_recursive(
    matrix<Type>& H,
    int& col_idx,
    const vector<Type>& current_col,
    int var_idx,
    int d,
    int remaining_deg,
    const matrix<Type>& InputMatrix
) {
    if (col_idx >= H.cols()) return;
    H.col(col_idx) = current_col;
    col_idx++;
    if (remaining_deg <= 0) return;
    for (int i = var_idx; i < d; i++) {
        vector<Type> next_col = (current_col.array() * InputMatrix.col(i).array()).matrix();
        fill_poly_cols_recursive(H, col_idx, next_col, i, d, remaining_deg - 1, InputMatrix);
    }
}

template <class Type>
void fill_poly(matrix<Type>& H, const matrix<Type>& InputMatrix, int deg) {
    int col_pos = 0;
    int n = InputMatrix.rows();
    int d = InputMatrix.cols();
    vector<Type> start_col(n);
    start_col.setOnes();
    fill_poly_cols_recursive(H, col_pos, start_col, 0, d, deg, InputMatrix);
}

// Function name MUST match the filename (poly_model)
template<class Type>
Type poly_model(objective_function<Type>* obj) {
  DATA_MATRIX(X);
  DATA_MATRIX(Y);
  DATA_INTEGER(deg);
  DATA_VECTOR(weights);
  PARAMETER_MATRIX(Theta);

  int n = X.rows();
  int p = Theta.rows();
  int layers = weights.size();

  matrix<Type> CurrentInput = X;
  matrix<Type> H(n, p);

  Type nll = 0;
  for(int k = 0; k < layers; k++) {
    fill_poly(H, CurrentInput, deg);
    CurrentInput = H * Theta;
    matrix<Type> diff = CurrentInput.topRows(n-k-1) - Y.bottomRows(n-k-1);
    nll += weights[k] * (diff.array() * diff.array()).sum();
  }
  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this