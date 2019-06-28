
#ifndef SVD_SOLVER_TRAITS_HPP_
#define SVD_SOLVER_TRAITS_HPP_

#include "svd_fwd.hpp"
#include "../../algebra/src/matrix/algebra_matrix_meta.hpp"
#include "../../algebra/src/multi_vector/algebra_multi_vector_meta.hpp"

namespace rompp{ 
namespace svd{
namespace details{

template<typename T, typename enable = void>
struct svd_traits{
  using derived_t = void;
  using matrix_t = void;
  using native_matrix_t = void;
  using scalar_t = void;
  using lsv_t = void;
  using rsv_t = void;
  using sval_t = void;
  
};
//---------------------------------------------------------------

template<typename T> 
struct svd_traits<const T> : svd_traits<T> {};
//---------------------------------------------------------------

  
#ifdef HAVE_TRILINOS
template <typename matrix_type,
	  template <typename...> class lsv_type,
	  template <typename...> class rsv_type,
	  typename sval_type>
struct svd_traits<Solver<
		    matrix_type,
		    lsv_type,
		    rsv_type,
		    sval_type,
		    typename
		    std::enable_if<
		      algebra::meta::is_sparse_matrix_epetra<
			typename
			algebra::details::traits<matrix_type>::wrapped_t
			>::value
		      >::type
		    >
		  >{

  using derived_t = Solver<matrix_type, lsv_type, rsv_type, sval_type>;

  using matrix_t = matrix_type;
  using native_matrix_t =
    typename algebra::details::traits<matrix_type>::wrapped_t;
  using scalar_t =
    typename algebra::details::traits<matrix_type>::scalar_t;
  using lsv_t = lsv_type<Epetra_MultiVector>;
  using rsv_t = rsv_type<Epetra_MultiVector>;
  using sval_t = sval_type;

};
#endif
//---------------------------------------------------------------

  
#ifdef HAVE_TRILINOS
template <typename matrix_type,
	  template <typename...> class lsv_type,
	  template <typename...> class rsv_type,
	  typename sval_type>
struct svd_traits<Solver<
		    matrix_type,
		    lsv_type,
		    rsv_type,
		    sval_type,
		    typename
		    std::enable_if<
		      algebra::meta::is_multi_vector_epetra<
			typename algebra::details::traits<matrix_type>::wrapped_t
			>::value
		      >::type
		    >
		  >{

  using derived_t = Solver<matrix_type, lsv_type, rsv_type, sval_type>;

  using matrix_t = matrix_type;
  using native_matrix_t =
    typename algebra::details::traits<matrix_type>::wrapped_t;
  using scalar_t =
    typename algebra::details::traits<matrix_type>::scalar_t;
  using lsv_t = lsv_type<Epetra_MultiVector>;
  using rsv_t = rsv_type<Epetra_MultiVector>;
  using sval_t = sval_type;
};
#endif

  
}//end namespace details
}//end namespace svd 
}//end namespace rompp
#endif
