#include "fast_mass_springs_precomputation_dense.h"
#include "signed_incidence_matrix_dense.h"
#include <Eigen/Dense>
#include <cmath>

bool fast_mass_springs_precomputation_dense(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & E,
  const double k,
  const Eigen::VectorXd & m,
  const Eigen::VectorXi & b,
  const double delta_t,
  Eigen::VectorXd & r,
  Eigen::MatrixXd & M,
  Eigen::MatrixXd & A,
  Eigen::MatrixXd & C,
  Eigen::LLT<Eigen::MatrixXd> & prefactorization)
{
  /////////////////////////////////////////////////////////////////////////////
  // Replace with your code
  Eigen::MatrixXd Q = Eigen::MatrixXd::Identity(V.rows(),V.rows());

  r.resize(E.rows());
  for (int idx = 0; idx < E.rows(); idx++){
    r(idx) = (V.row(E(idx, 0)) - V.row(E(idx, 1))).norm();
  }

  // the diagonal matrix M
  M = Eigen::MatrixXd::Zero(V.rows(), V.rows());
  for (int idx = 0; idx < M.rows(); idx++){
		M(idx, idx) = m(idx);
  }

  // calculate matrix A
  A = Eigen::MatrixXd::Zero(E.rows(), V.rows());
  signed_incidence_matrix_dense(V.rows(), E, A);

  // calculate C
  C = Eigen::MatrixXd::Zero(b.rows(), V.rows()); 
  for (int idx = 0; idx < b.rows(); idx++){
    C(idx, b(idx)) = 1;
  }

  // use the formula Q = k*A^T*A + 1/(delta_t^2) * M
  Q = k * A.transpose() * A + 1/ (pow(delta_t, 2)) * M;

	double w = 1e10;  
  Q	+= w * C.transpose() * C;
  /////////////////////////////////////////////////////////////////////////////
  prefactorization.compute(Q);
  return prefactorization.info() != Eigen::NumericalIssue;
}
