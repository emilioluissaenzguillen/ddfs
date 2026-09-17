/* Routines for ddfs R-package
 *
 * Routines are based on the algorithm in Dimitrova et al. (2025).
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, a copy is available at
 * https://www.R-project.org/Licenses/
 *
 * These functions are distributed WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.	 See the GNU General Public License for more details.
 */

#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector constrMLE(NumericVector theta,
                        const NumericVector indexes,
                        const NumericMatrix basisMatrix,
                        const NumericVector first_term) {

  NumericVector theta_prev = clone(theta);
  NumericVector f_X_hat_prev(basisMatrix.nrow());
  bool ok = true;

  while (ok) {
    // Update theta_prev only for the given indexes
    for (int i = 0; i < indexes.size(); i++) {
      int idx = indexes[i] - 1; // Adjust for 0-based indexing in C++
      theta_prev[idx] = theta[idx];
    }

    // Compute f_X_hat_prev = basisMatrix %*% theta_prev
    for (int i = 0; i < basisMatrix.nrow(); i++) {
      // Extract the i-th row as a NumericVector
      NumericVector row = basisMatrix.row(i);
      // Compute the inner product of this row and theta_prev
      f_X_hat_prev[i] = std::inner_product(row.begin(), row.end(), theta_prev.begin(), 0.0);
    }

    // Ensure no zero values in f_X_hat_prev
    for (int i = 0; i < f_X_hat_prev.size(); i++) {
      if (f_X_hat_prev[i] == 0) f_X_hat_prev[i] = 1e-6;
    }

    // Compute theta update: theta = first_term * theta_prev * colSums( basisMatrix / f_X_hat_prev )
    NumericVector col_sums(basisMatrix.ncol(), 0.0);

    for (int j = 0; j < basisMatrix.ncol(); j++) {
      for (int i = 0; i < basisMatrix.nrow(); i++) {
        col_sums[j] += basisMatrix(i, j) / f_X_hat_prev[i];
      }
    }

    for (int i = 0; i < theta.size(); i++) {
      theta[i] = first_term[i] * theta_prev[i] * col_sums[i];
    }

    // Convergence check
    double max_diff = 0.0;
    for (int i = 0; i < indexes.size(); i++) {
      int idx = indexes[i] - 1; // Adjust for 0-based indexing
      double diff = std::abs(theta[idx] - theta_prev[idx]);
      if (diff > max_diff) max_diff = diff;
    }

    if (max_diff < 1e-6) ok = false;
  }

  return theta;
}


// [[Rcpp::export]]
NumericVector constrMLE_ef(NumericVector theta,
                           const NumericVector indexes,
                           const NumericMatrix basisMatrix,
                           const NumericVector first_term) {
  int nrow = basisMatrix.nrow();
  int ncol = basisMatrix.ncol();
  int ntheta = theta.size();

  // Precompute adjusted indexes (convert from 1-indexed to 0-indexed)
  std::vector<int> idx_vec(indexes.size());
  for (int i = 0; i < indexes.size(); i++){
    idx_vec[i] = indexes[i] - 1;
  }

  NumericVector theta_prev = clone(theta);
  NumericVector f_X_hat_prev(nrow);
  bool ok = true;

  while(ok) {
    // Update theta_prev for the given indexes
    for (size_t k = 0; k < idx_vec.size(); k++) {
      int idx = idx_vec[k];
      theta_prev[idx] = theta[idx];
    }

    // Compute f_X_hat_prev = basisMatrix %*% theta_prev without creating temporary vectors
    for (int i = 0; i < nrow; i++) {
      double sum = 0.0;
      for (int j = 0; j < ncol; j++){
        sum += basisMatrix(i, j) * theta_prev[j];
      }
      // Replace zeros immediately
      f_X_hat_prev[i] = (sum == 0.0) ? 1e-6 : sum;
    }

    // Compute col_sums for each column of basisMatrix
    NumericVector col_sums(ncol, 0.0);
    for (int j = 0; j < ncol; j++){
      double sum = 0.0;
      // Since basisMatrix is column-major, this loop is cache friendly
      for (int i = 0; i < nrow; i++){
        sum += basisMatrix(i, j) / f_X_hat_prev[i];
      }
      col_sums[j] = sum;
    }

    // Update theta
    for (int i = 0; i < ntheta; i++){
      theta[i] = first_term[i] * theta_prev[i] * col_sums[i];
    }

    // Convergence check: only consider the adjusted indexes
    double max_diff = 0.0;
    for (size_t k = 0; k < idx_vec.size(); k++){
      int idx = idx_vec[k];
      double diff = std::abs(theta[idx] - theta_prev[idx]);
      if (diff > max_diff) max_diff = diff;
    }

    if (max_diff < 1e-6) {
      ok = false;
    }
  }

  return theta;
}


// [[Rcpp::export]]
arma::rowvec constrMLE_arma(arma::vec theta,
                            arma::uvec indexes,
                            arma::mat basisMatrix,
                            arma::vec first_term) {
  // Adjust indexes from 1-indexed (R) to 0-indexed (C++)
  arma::uvec idx = indexes;
  idx -= 1;

  arma::vec theta_prev = theta;
  arma::vec f_X_hat_prev;
  bool ok = true;

  while(ok) {
    // Update theta_prev for the specified indexes
    theta_prev.elem(idx) = theta.elem(idx);

    // Compute f_X_hat_prev = basisMatrix %*% theta_prev
    f_X_hat_prev = basisMatrix * theta_prev;

    // Replace zeros in f_X_hat_prev with 1e-6
    arma::uvec zeroIdx = arma::find(f_X_hat_prev == 0);
    if(zeroIdx.n_elem > 0) {
      f_X_hat_prev.elem(zeroIdx).fill(1e-6);
    }

    // Compute column sums: equivalent to colSums(basisMatrix / f_X_hat_prev)
    arma::mat ratio = basisMatrix;
    ratio.each_col() /= f_X_hat_prev;
    arma::vec col_sums = arma::sum(ratio, 0).t(); // get column sums as column vector

    // Update theta element-wise: theta = first_term * theta_prev * col_sums
    theta = first_term % theta_prev % col_sums;

    // Convergence check for specified indices
    if (arma::max(arma::abs(theta.elem(idx) - theta_prev.elem(idx))) < 1e-6)
      ok = false;
  }

  // Return the result as a row vector (transpose the column vector theta)
  return arma::trans(theta);
}


