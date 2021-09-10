
#ifndef PRESSIO_TESTS_WLS_ROM_DATA_TYPE_KOKKOS_HPP_
#define PRESSIO_TESTS_WLS_ROM_DATA_TYPE_KOKKOS_HPP_

namespace {

template <typename scalar_t>
struct romDataTypeKokkos
{
  using execution_space = Kokkos::DefaultExecutionSpace;
  using kll		= Kokkos::LayoutLeft;
  using k1dLl_d		= Kokkos::View<scalar_t*, kll, execution_space>;
  using k2dLl_d		= Kokkos::View<scalar_t**, kll, execution_space>;
  using wls_state_t	= pressio::containers::Vector<k1dLl_d>;
  using wls_hessian_t	= pressio::containers::DenseMatrix<k2dLl_d>;
  using lin_solver_tag  = pressio::solvers::linear::direct::potrsL;
  using linear_solver_t = pressio::solvers::linear::Solver<lin_solver_tag, wls_hessian_t>;
};

}// end namespace

#endif
