/*
//@HEADER
// ************************************************************************
//
// qr_rfactor_solve_impl.hpp
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

#ifndef QR_IMPL_QR_RFACTOR_SOLVE_IMPL_HPP_
#define QR_IMPL_QR_RFACTOR_SOLVE_IMPL_HPP_

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseSolver.hpp"
#include <Teuchos_SerialQRDenseSolver.hpp>
#endif

namespace pressio{ namespace qr{ namespace impl{

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template< typename vector_type, typename R_type>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_vector_wrapper_eigen<vector_type>::value and
  ::pressio::containers::predicates::is_dense_matrix_wrapper_eigen<R_type>::value
>
solve(const vector_type & rhs,
	   const R_type & Rmatrix,
	   vector_type & y){
  auto nat_rhs = rhs.data();
  auto nat_y = y.data();
  auto nat_R = Rmatrix.data();

  *nat_y = (*nat_R).template triangularView<Eigen::Upper>().solve(*nat_rhs);
}
#endif

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
template<typename vector_type, typename R_type>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_vector_wrapper_teuchos<vector_type>::value and
  ::pressio::containers::predicates::is_dense_matrix_teuchos_rcp<R_type>::value
>
solve(const vector_type & rhs, R_type Rmatrix, vector_type & y)
{
  using ord_t	 = typename R_type::element_type::ordinalType;
  using sc_t	 = typename R_type::element_type::scalarType;
  using wrapv_t	 = typename containers::details::traits<vector_type>::wrapped_t;
  using solver_t = Teuchos::SerialDenseSolver<ord_t, sc_t>;

  solver_t My_Solver;
  My_Solver.setMatrix(Rmatrix);
  My_Solver.setVectors( Teuchos::rcp( y.data(), false ),
			Teuchos::rcp( const_cast<wrapv_t*>(rhs.data()), false ) );

  int info = My_Solver.factor();
  if (info != 0)
    std::cout << "SerialDenseSolver::factor() returned : "
	      << info << std::endl;

  info = My_Solver.solve();
  if (info != 0)
    std::cout << "SerialDenseSolver::solve() returned : "
	      << info << std::endl;
}
#endif

#if defined(PRESSIO_ENABLE_TPL_TRILINOS)
template<typename vector_type, typename R_type>
::pressio::mpl::enable_if_t<
  !::pressio::containers::predicates::is_vector_wrapper_teuchos<vector_type>::value and
  ::pressio::containers::predicates::is_dense_matrix_teuchos_rcp<R_type>::value
>
solve(const vector_type & rhs, R_type Rmatrix, vector_type & y)
{
  // todo: maybe here we should use directly backward substitution but it does
  // not seem to be available directly from teuchos

  using ord_t	  = typename R_type::element_type::ordinalType;
  using sc_t	  = typename R_type::element_type::scalarType;
  using tservec_t = Teuchos::SerialDenseVector<ord_t, sc_t>;

  auto vecSize = rhs.extent(0);
  tservec_t rhsTV(Teuchos::View, const_cast<sc_t*>(rhs.data()->data()), vecSize);
  tservec_t yTV(Teuchos::View, y.data()->data(), vecSize);

  Teuchos::SerialQRDenseSolver<ord_t, sc_t> My_Solver;
  My_Solver.setMatrix(Rmatrix);
  My_Solver.setVectors( Teuchos::rcp(&yTV, false), Teuchos::rcp(&rhsTV, false) );
  int info = My_Solver.factor();
  if (info != 0)
    std::cout << "SerialDenseSolver::factor() returned : " << info << std::endl;

  info = My_Solver.solve();
  if (info != 0)
    std::cout << "SerialDenseSolver::solve() returned : " << info << std::endl;
}
#endif

}}}//end namespace pressio::qr::impl
#endif  // QR_IMPL_QR_RFACTOR_SOLVE_IMPL_HPP_
