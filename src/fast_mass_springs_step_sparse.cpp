#include "fast_mass_springs_step_sparse.h"
#include <igl/matlab_format.h>

void fast_mass_springs_step_sparse(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & E,
  const double k,
  const Eigen::VectorXi & b,
  const double delta_t,
  const Eigen::MatrixXd & fext,
  const Eigen::VectorXd & r,
  const Eigen::SparseMatrix<double>  & M,
  const Eigen::SparseMatrix<double>  & A,
  const Eigen::SparseMatrix<double>  & C,
  const Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > & prefactorization,
  const Eigen::MatrixXd & Uprev,
  const Eigen::MatrixXd & Ucur,
  Eigen::MatrixXd & Unext)
{
  //////////////////////////////////////////////////////////////////////////////
  // Replace with your code
  Eigen::MatrixXd d = Eigen::MatrixXd::Zero(E.rows(), 3);
  Eigen::MatrixXd p = Ucur;

  for(int iter = 0;iter < 50;iter++)
  {
    for (int idx=0; idx<d.rows(); idx++){
			d.row(idx) = r(idx) * (p.row(E(idx,0)) - p.row(E(idx, 1))).normalized();
 	 	}
      
    Eigen::MatrixXd y = 1/ (pow(delta_t, 2)) * M * (2 * Ucur - Uprev) + fext;
    const Eigen::MatrixXd l = k * A.transpose() * d + y + 1e10 * C.transpose() * C * V;
    p = prefactorization.solve(l);
  }
  Unext = p;
  //////////////////////////////////////////////////////////////////////////////
}
