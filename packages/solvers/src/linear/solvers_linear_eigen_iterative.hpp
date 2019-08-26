
#ifndef SOLVERS_LINEAR_EIGEN_ITERATIVE_HPP
#define SOLVERS_LINEAR_EIGEN_ITERATIVE_HPP

#include "../solvers_ConfigDefs.hpp"
#include "../base/solvers_iterative_base.hpp"
#include "../base/solvers_linear_base.hpp"
#include "solvers_linear_traits.hpp"

namespace pressio { namespace solvers { namespace iterative{

template<typename SolverT, typename MatrixT>
class EigenIterative
  : public LinearBase<SolverT, MatrixT, EigenIterative<SolverT, MatrixT> >,
    public IterativeBase< EigenIterative<SolverT, MatrixT>,
			  typename containers::details::traits<MatrixT>::scalar_t>
{
public:

  static_assert( ::pressio::containers::meta::is_matrix_wrapper_eigen<MatrixT>::value or
  		 ::pressio::containers::meta::is_multi_vector_wrapper_eigen<MatrixT>::value,
  		 "Eigen iterative solver needs a matrix type = wrapper of an eigen matrix");

  using solver_t	= SolverT;
  using matrix_type	= MatrixT;
  using native_mat_t    = typename containers::details::traits<MatrixT>::wrapped_t;
  using scalar_t        = typename containers::details::traits<MatrixT>::scalar_t;
  using this_t          = EigenIterative<SolverT, MatrixT>;
  using base_interface  = LinearBase<SolverT, MatrixT, this_t>;
  using base_iterative  = IterativeBase<this_t, scalar_t>;
  using solver_traits   = linear::details::traits<SolverT>;
  using native_solver_t = typename solver_traits::template eigen_solver_type<native_mat_t>;
  using iteration_t	= typename base_iterative::iteration_t;

  static_assert( solver_traits::eigen_enabled == true,
		 "the native solver must be from Eigen to use in EigenIterative");

  static_assert( solver_traits::direct == false,
		 "The native eigen solver must be iterative to use in EigenIterative");

public:
  EigenIterative() = default;
  EigenIterative(const EigenIterative &) = delete;
  ~EigenIterative() = default;

private:

  iteration_t getNumIterationsExecutedImpl() const {
    return mysolver_.iterations();
  }

  scalar_t getFinalErrorImpl() const {
    return mysolver_.error();
  }

  void resetLinearSystemImpl(const MatrixT& A) {
    mysolver_.compute(*A.data());
  }

  template <typename T>
  void solveImpl(const T& b, T & y)
  {
    *y.data() = mysolver_.solve(*b.data());
  }

  template <typename T>
  void solveImpl(const MatrixT & A, const T& b, T & y) {
    this->resetLinearSystem(A);
    this->solve(b, y);
  }

  template <typename T>
  void solveAllowMatOverwriteImpl(MatrixT & A, const T& b, T & y) {
    this->resetLinearSystem(A);
    this->solve(b, y);
  }

  friend base_iterative;
  friend base_interface;
  native_solver_t mysolver_ = {};
};

}}} // end namespace pressio::solvers::iterative
#endif
