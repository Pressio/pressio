#ifndef ROM_DATA_TYPE_EIGEN_HPP_
#define ROM_DATA_TYPE_EIGEN_HPP_

template <typename scalar_t>
struct romDataTypeEigen
{
  using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
  using eig_dyn_mat	= Eigen::Matrix<scalar_t, -1, -1>;
  using wls_state_t	= pressio::containers::Vector<eig_dyn_vec>;
  using wls_hessian_t   = pressio::containers::Matrix<eig_dyn_mat>;
  using lin_solver_tag	= pressio::solvers::linear::direct::potrsL;
  using linear_solver_t = pressio::solvers::linear::Solver<lin_solver_tag, wls_hessian_t>;
};

// end namespace
#endif
