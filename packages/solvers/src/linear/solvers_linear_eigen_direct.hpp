
#ifndef SOLVERS_LINEAR_EIGEN_DIRECT_HPP
#define SOLVERS_LINEAR_EIGEN_DIRECT_HPP

#include "../solvers_ConfigDefs.hpp"
#include "../base/solvers_linear_base.hpp"
#include "solvers_linear_traits.hpp"

namespace rompp { namespace solvers { namespace direct{

template<typename SolverT, typename MatrixT>
class EigenDirect
  : public LinearBase<SolverT, MatrixT,
		      EigenDirect<SolverT, MatrixT> >
{
public:
  static_assert( ::rompp::core::meta::is_matrix_wrapper_eigen<MatrixT>::value or
  		 ::rompp::core::meta::is_multi_vector_wrapper_eigen<MatrixT>::value,
  		 "Eigen direct solver needs a matrix type = wrapper of an eigen matrix");

  using solver_t	= SolverT;
  using matrix_type	= MatrixT;
  using native_mat_t    = typename core::details::traits<MatrixT>::wrapped_t;
  using scalar_t        = typename core::details::traits<MatrixT>::scalar_t;
  using this_t          = EigenDirect<SolverT, MatrixT>;
  using base_interface  = LinearBase<SolverT, MatrixT, this_t>;
  using solver_traits   = linear::details::traits<SolverT>;
  using native_solver_t = typename solver_traits::template eigen_solver_type<native_mat_t>;

  static_assert( solver_traits::eigen_enabled == true,
		 "the native solver must be from Eigen to use in EigenDirect");

  static_assert( solver_traits::direct == true,
		 "the native eigen solver must be direct to use in EigenDirect");

public:
  EigenDirect() = default;
  EigenDirect(const EigenDirect &) = delete;
  ~EigenDirect() = default;

private:
  void resetLinearSystemImpl(const MatrixT& A) {
    mysolver_.compute(*A.data());
  }

  template <typename T>
  void solveImpl(const T& b, T & y) {
    *y.data() = mysolver_.solve(*b.data());
  }

  friend base_interface;
  native_solver_t mysolver_ = {};
};

}}} // end namespace rompp::solvers::direct
#endif
