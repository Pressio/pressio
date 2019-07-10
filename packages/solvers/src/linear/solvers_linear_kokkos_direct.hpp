#ifdef HAVE_TRILINOS
#ifndef SOLVERS_LINEAR_KOKKOS_DIRECT_HPP
#define SOLVERS_LINEAR_KOKKOS_DIRECT_HPP

#include "../solvers_ConfigDefs.hpp"
#include "../base/solvers_linear_base.hpp"
#include "solvers_linear_traits.hpp"
#include <Teuchos_LAPACK.hpp>
#include <Teuchos_SerialDenseSolver.hpp>

namespace pressio { namespace solvers { namespace direct{

template<typename SolverT, typename MatrixT, typename enable = void>
class KokkosDirect;

/* GETRS */
template<typename SolverT, typename MatrixT>
class KokkosDirect<
  SolverT, MatrixT,
  mpl::enable_if_t<
    mpl::is_same<SolverT, ::pressio::solvers::linear::direct::getrs>::value
    >
  >
  : public LinearBase<SolverT, MatrixT, KokkosDirect<SolverT, MatrixT>>
{
public:
  static_assert( ::pressio::containers::meta::is_dense_matrix_wrapper_kokkos<MatrixT>::value or
  		 ::pressio::containers::meta::is_multi_vector_wrapper_kokkos<MatrixT>::value,
  		 "Kokkos direct dense solver expects either (a) dense matrix wrapper or a (b) multi-vector wrapper, both wrapping a rank=2 Kokkos View");

  using this_t          = KokkosDirect<SolverT, MatrixT>;
  using solver_t	= SolverT;
  using matrix_type	= MatrixT;
  using native_mat_t    = typename containers::details::traits<MatrixT>::wrapped_t;
  using scalar_t        = typename containers::details::traits<MatrixT>::scalar_t;
  using exe_space       = typename containers::details::traits<MatrixT>::execution_space;

  using base_t  = LinearBase<SolverT, MatrixT, this_t>;
  using solver_traits   = linear::details::traits<SolverT>;

  static_assert( solver_traits::kokkos_enabled == true,
  		 "the native solver must suppport kokkos to use in KokkosDirect");
  static_assert( solver_traits::direct == true,
  		 "the native solver must be direct to use in KokkosDirect");

public:
  KokkosDirect() = default;
  KokkosDirect(const KokkosDirect &) = delete;
  ~KokkosDirect() = default;

private:

  /*
   * enable if:
   * the matrix has layout left (i.e. column major)
   * T is a kokkos vector wrapper
   * has host execution space
   * T and MatrixT have same execution space
   */
  template <
    typename _MatrixT = MatrixT,
    typename _SolverT = SolverT,
    typename T,
    mpl::enable_if_t<
      mpl::is_same<
	typename ::pressio::containers::details::traits<_MatrixT>::layout,
	Kokkos::LayoutLeft
	>::value
      and
      ::pressio::containers::meta::is_vector_wrapper_kokkos<T>::value
      and
#ifdef KOKKOS_ENABLE_OPENMP
      (
#endif
       mpl::is_same<
       typename containers::details::traits<T>::execution_space,
       Kokkos::Serial
       >::value
#ifdef KOKKOS_ENABLE_OPENMP
       or
       mpl::is_same<
       typename containers::details::traits<T>::execution_space,
       Kokkos::OpenMP
       >::value)
#endif
	and
	mpl::is_same<
	typename containers::details::traits<T>::execution_space,
      typename containers::details::traits<_MatrixT>::execution_space
      >::value
      > * = nullptr
  >
  void solveAllowMatOverwriteImpl(_MatrixT & A, const T& b, T & y) {
    assert(A.rows() == b.size() );
    assert(A.cols() == y.size() );

    // gerts is for square matrices
    assert(A.rows() == A.cols() );

    const auto nRows = A.rows();
    const auto nCols = A.cols();

    int info = 0;
    const int ipivSz = std::min(nRows, nCols);
    int ipiv[ipivSz];

    // LU factorize using GETRF
    lpk_.GETRF(nRows, nCols, A.data()->data(),
	      nRows, // leading dim is nRows because it is col-major
	      ipiv,
	      &info);
    assert(info == 0);

    Kokkos::deep_copy(*y.data(), *b.data());

    const char ct = 'N';
    const int nrhs = 1; // only solving for a single vector b
    lpk_.GETRS(ct,
	      nCols,
    	      nrhs,
    	      A.data()->data(),
    	      nCols,
	      ipiv,
    	      y.data()->data(),
    	      y.size(),
    	      &info);
    assert(info == 0);
  }

  friend base_t;
  Teuchos::LAPACK<int, scalar_t> lpk_;
};

}}} // end namespace pressio::solvers::direct
#endif
#endif
