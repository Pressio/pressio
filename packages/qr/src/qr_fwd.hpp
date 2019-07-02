
#ifndef QR_FORWARD_DECLARATIONS_HPP_
#define QR_FORWARD_DECLARATIONS_HPP_

#include "qr_ConfigDefs.hpp"
#include "../../containers/src/containers_fwd.hpp"

namespace pressio{  namespace qr{

// all fwd decl of mixin base classes
template<typename derived, typename matrix_t>
class QRInPlaceBase;

template<typename derived,
	 typename matrix_t, typename Q_type>
class QROutOfPlaceBase;

template<typename derived, typename R_type>
class RFactorBase;

template<typename derived>
class QRSolveBase;


namespace impl{

template<
  typename matrix_type,
  typename R_t = void,
  int n = utils::constants::dynamic,
  int m = utils::constants::dynamic,
  template <typename...> class Q_type
	= containers::MultiVector,
  typename enable= void>
class QRHouseholderDenseEigenMatrixWrapper;

#if defined HAVE_TRILINOS
template<typename matrix_t,
	 typename R_t,
	 int n = utils::constants::dynamic,
	 int m = utils::constants::dynamic,
	 template <typename...> class Q_type
	 = containers::MultiVector>
class EpetraMVHouseholderUsingEigen;

template<typename matrix_t,
	 typename R_t,
	 int n = utils::constants::dynamic,
	 int m = utils::constants::dynamic,
	 template <typename...> class Q_type
	 = containers::MultiVector>
class TpetraMVHouseholderUsingEigen;

template<typename matrix_t,
	 typename R_t,
	 int n, int m,
	 typename wrap_Q_type,
	 template <typename...> class Q_type
	 = containers::MultiVector,
	 typename enable = void>
class EpetraMVTSQR;

template<typename matrix_t,
	 typename R_t,
	 int n, int m,
	 typename wrap_Q_type,
	 template <typename...> class Q_type
	 = containers::MultiVector,
	 typename enable = void>
class TpetraMVTSQR;

template<typename matrix_t,
	 typename R_t,
	 int n, int m,
	 typename wrap_Q_type,
	 template <typename...> class Q_type
	 = containers::MultiVector,
	 typename enable = void>
class TpetraBlockMVTSQR;

template<typename matrix_t,
	 typename R_t,
	 int n = utils::constants::dynamic,
	 int m = utils::constants::dynamic,
	 typename wrap_Q_type = void,
	 template <typename...> class Q_type
	 = containers::MultiVector,
	 typename enable = void>
class ModGramSchmidtMVEpetra;

template<typename matrix_t,
	 typename R_t,
	 int n = utils::constants::dynamic,
	 int m = utils::constants::dynamic,
	 typename wrap_Q_type = void,
	 template <typename...> class Q_type
	 = containers::MultiVector,
	 typename enable = void>
class ModGramSchmidtMVTpetra;

#endif //HAVE_TRILINOS


template<
  typename matrix_type,
  typename algorithm,
  bool in_place,
  int m,
  int n,
  typename R_type,
  template <typename...> class Q_type,
  typename enable = void>
class QRSolver;

}//end namespace pressio::qr::impl


template<
  typename matrix_type,
  typename algorithm,
  bool in_place = false,
  int n = utils::constants::dynamic,
  int m = utils::constants::dynamic,
  template <typename...> class Q_type
        = ::pressio::containers::MultiVector,
  typename enable = void
  >
using QRSolver = impl::QRSolver<matrix_type, algorithm,
				in_place, m, n,
				void, Q_type>;

template<
  typename matrix_type,
  typename algorithm,
  typename R_type,
  bool in_place = false,
  int n = utils::constants::dynamic,
  int m = utils::constants::dynamic,
  template <typename...> class Q_type
        = ::pressio::containers::MultiVector,
  typename enable = void
  >
using QRSolverWrapR = impl::QRSolver<matrix_type,
				     algorithm,
				     in_place,
				     m, n,
				     R_type,
				     Q_type>;


}}//end namespace pressio::qr
#endif
