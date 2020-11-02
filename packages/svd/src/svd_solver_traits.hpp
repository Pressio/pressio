/*
//@HEADER
// ************************************************************************
//
// svd_solver_traits.hpp
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef SVD_SVD_SOLVER_TRAITS_HPP_
#define SVD_SVD_SOLVER_TRAITS_HPP_

namespace pressio{  namespace svd{ namespace details{

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

  
// #ifdef PRESSIO_ENABLE_TPL_TRILINOS
// template <typename matrix_type,
// 	  template <typename...> class lsv_type,
// 	  template <typename...> class rsv_type,
// 	  typename sval_type>
// struct svd_traits<Solver<
// 		    matrix_type,
// 		    lsv_type,
// 		    rsv_type,
// 		    sval_type,
// 		    typename
// 		    std::enable_if<
// 		      containers::predicates::is_sparse_matrix_epetra<
// 			typename
// 			containers::details::traits<matrix_type>::wrapped_t
// 			>::value
// 		      >::type
// 		    >
// 		  >{

//   using derived_t = Solver<matrix_type, lsv_type, rsv_type, sval_type>;

//   using matrix_t = matrix_type;
//   using native_matrix_t =
//     typename containers::details::traits<matrix_type>::wrapped_t;
//   using scalar_t =
//     typename containers::details::traits<matrix_type>::scalar_t;
//   using lsv_t = lsv_type<Epetra_MultiVector>;
//   using rsv_t = rsv_type<Epetra_MultiVector>;
//   using sval_t = sval_type;

// };
// #endif
// //---------------------------------------------------------------

  
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
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
		      containers::predicates::is_multi_vector_epetra<
			typename containers::details::traits<matrix_type>::wrapped_t
			>::value
		      >::type
		    >
		  >{

  using derived_t = Solver<matrix_type, lsv_type, rsv_type, sval_type>;

  using matrix_t = matrix_type;
  using native_matrix_t =
    typename containers::details::traits<matrix_type>::wrapped_t;
  using scalar_t =
    typename containers::details::traits<matrix_type>::scalar_t;
  using lsv_t = lsv_type<Epetra_MultiVector>;
  using rsv_t = rsv_type<Epetra_MultiVector>;
  using sval_t = sval_type;
};
#endif

  
}//end namespace details
}//end namespace svd 
}//end namespace pressio
#endif  // SVD_SVD_SOLVER_TRAITS_HPP_
