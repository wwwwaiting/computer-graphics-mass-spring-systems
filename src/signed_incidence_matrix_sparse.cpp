#include "signed_incidence_matrix_sparse.h"
#include <vector>

void signed_incidence_matrix_sparse(
  const int n,
  const Eigen::MatrixXi & E,
  Eigen::SparseMatrix<double>  & A)
{
  //////////////////////////////////////////////////////////////////////////////
  // Replace with your code
  std::vector<Eigen::Triplet<double> > ijv;
  for (int idx = 0; idx < E.rows(); idx++){
    ijv.emplace_back(Eigen::Triplet<double>(idx, E(idx, 0), 1));
    ijv.emplace_back(Eigen::Triplet<double>(idx, E(idx, 1), -1));
  }
  A.resize(E.rows(),n);
  A.setFromTriplets(ijv.begin(),ijv.end());
  //////////////////////////////////////////////////////////////////////////////
}
