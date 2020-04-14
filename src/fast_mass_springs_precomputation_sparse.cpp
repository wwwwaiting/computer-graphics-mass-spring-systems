#include "fast_mass_springs_precomputation_sparse.h"
#include "signed_incidence_matrix_sparse.h"
#include <vector>
#include <cmath>

bool fast_mass_springs_precomputation_sparse(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & E,
  const double k,
  const Eigen::VectorXd & m,
  const Eigen::VectorXi & b,
  const double delta_t,
  Eigen::VectorXd & r,
  Eigen::SparseMatrix<double>  & M,
  Eigen::SparseMatrix<double>  & A,
  Eigen::SparseMatrix<double>  & C,
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > & prefactorization)
{
  /////////////////////////////////////////////////////////////////////////////
  // Replace with your code
  const int n = V.rows();
  Eigen::SparseMatrix<double> Q(n,n);

  r.resize(E.rows());
  for (int idx = 0; idx < E.rows(); idx++){
    r(idx) = (V.row(E(idx, 0)) - V.row(E(idx, 1))).norm();
  }

 // the diagonal matrix M
  M.resize(n, n);
  std::vector<Eigen::Triplet<double> > ijv_m;
  for (int idx=0; idx<M.rows(); idx++){
  	ijv_m.emplace_back(idx,idx,m(idx));
  }
  M.setFromTriplets(ijv_m.begin(),ijv_m.end());

  // get matrix A
  signed_incidence_matrix_sparse(V.rows(), E, A);

  // calculate C
  C.resize(b.rows(), n);
  std::vector<Eigen::Triplet<double> > ijv_c;
  for (int idx=0; idx<C.rows(); idx++){
  	ijv_c.emplace_back(idx,b(idx),1);
  }
  C.setFromTriplets(ijv_c.begin(),ijv_c.end());

  // use the formula Q = k*A^T*A + 1/(delta_t^2) * M
  Q = k * A.transpose() * A + 1/ (pow(delta_t, 2)) * M;

	double w = 1e10;  
  Q	+= w * C.transpose() * C;

  /////////////////////////////////////////////////////////////////////////////
  prefactorization.compute(Q);
  return prefactorization.info() != Eigen::NumericalIssue;
}
