
#ifndef SOLVERS_LINEAR_SOLVER_IMPL_HPP_
#define SOLVERS_LINEAR_SOLVER_IMPL_HPP_

#include "solvers_linear_eigen_direct_impl.hpp"
#include "solvers_linear_eigen_iterative_impl.hpp"
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
#include "solvers_linear_kokkos_direct_geqrf_impl.hpp"
#include "solvers_linear_kokkos_direct_getrs_impl.hpp"
#include "solvers_linear_kokkos_direct_potrs_lower_impl.hpp"
#include "solvers_linear_kokkos_direct_potrs_upper_impl.hpp"
#endif

namespace pressio{ namespace solvers{ namespace linear { namespace impl{

template<typename tag, typename MatrixT, typename enable = void>
struct LinearSolverSelector{
  using type = void;
};

template<typename tag, typename MatrixT>
struct LinearSolverSelector<
  tag, MatrixT,
  mpl::enable_if_t<
    ::pressio::solvers::linear::details::traits<tag>::iterative and
    (::pressio::containers::meta::is_matrix_wrapper_eigen<MatrixT>::value or
     ::pressio::containers::meta::is_multi_vector_wrapper_eigen<MatrixT>::value)
    >
  >
{
  using native_mat_t    = typename containers::details::traits<MatrixT>::wrapped_t;
  using solver_traits   = linear::details::traits<tag>;
  using type = ::pressio::solvers::linear::impl::EigenIterative<tag, MatrixT>;
};

template<typename tag, typename MatrixT>
struct LinearSolverSelector<
  tag, MatrixT,
  mpl::enable_if_t<
    ::pressio::solvers::linear::details::traits<tag>::direct and
    (::pressio::containers::meta::is_matrix_wrapper_eigen<MatrixT>::value or
     ::pressio::containers::meta::is_multi_vector_wrapper_eigen<MatrixT>::value)
    >
  >
{
  using native_mat_t    = typename containers::details::traits<MatrixT>::wrapped_t;
  using solver_traits   = linear::details::traits<tag>;
  using type = ::pressio::solvers::linear::impl::EigenDirect<tag, MatrixT>;
};


#ifdef PRESSIO_ENABLE_TPL_KOKKOS
template<typename tag, typename MatrixT>
struct LinearSolverSelector<
  tag, MatrixT,
  mpl::enable_if_t<
    ::pressio::solvers::linear::details::traits<tag>::direct and
    (::pressio::containers::meta::is_dense_matrix_wrapper_kokkos<MatrixT>::value or
     ::pressio::containers::meta::is_multi_vector_wrapper_kokkos<MatrixT>::value)
    >
  >
{
  using native_mat_t    = typename containers::details::traits<MatrixT>::wrapped_t;
  using solver_traits   = linear::details::traits<tag>;
  using type = ::pressio::solvers::linear::impl::KokkosDirect<tag, MatrixT>;
};
#endif

}}}}// end namespace pressio::solvers::linear::impl
#endif
