/*
//@HEADER
// ************************************************************************
//
// ops_level3.hpp
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

#ifndef OPS_PYBIND11_OPS_LEVEL3_HPP_
#define OPS_PYBIND11_OPS_LEVEL3_HPP_

namespace pressio{ namespace ops{

/*
 * C = beta * C + alpha*A*B
 * A,B,C = rank-2 tensors
*/
template <typename A_type, typename B_type, typename scalar_type, typename C_type>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<A_type>::value and
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<B_type>::value and
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<C_type>::value
  >
product(::pressio::nontranspose modeA,
	::pressio::nontranspose modeB,
	const scalar_type alpha,
	const A_type & A,
	const B_type & B,
	const scalar_type beta,
	C_type & C)
{
  static_assert
    (containers::predicates::are_scalar_compatible<A_type, B_type, C_type>::value,
     "Types are not scalar compatible");

  // NOTE: need to check if doing this import is expensive,
  // and assess whether we can use blas directly when we know
  // that objects involved are dense with not strange layout.
  pybind11::object pyblas = pybind11::module::import("scipy.linalg.blas");
  constexpr auto one  = ::pressio::utils::constants<scalar_type>::one();
  constexpr auto no   = ::pressio::utils::constants<int>::zero();
  constexpr auto yes  = ::pressio::utils::constants<int>::one();
  constexpr auto transA = no;
  constexpr auto transB = no;
  constexpr auto ovw    = yes;
  pyblas.attr("dgemm")(one, *A.data(), *B.data(), beta,
		       *C.data(), transA, transB, ovw);
}

/*
 * C = beta * C + alpha*A*B
 * where A,B,C = rank-3 tensor
*/
template <
  typename A_type, typename B_type, typename scalar_type, typename C_type
  >
::pressio::mpl::enable_if_t<
  containers::predicates::is_fstyle_array_pybind<A_type>::value and
  containers::predicates::is_fstyle_array_pybind<B_type>::value and
  containers::predicates::is_fstyle_array_pybind<C_type>::value
  >
product(::pressio::nontranspose,
	::pressio::nontranspose,
	const scalar_type alpha,
	const ::pressio::containers::Tensor<3, A_type> & A,
	const ::pressio::containers::Tensor<3, B_type> & B,
	const scalar_type beta,
	::pressio::containers::Tensor<3, C_type> & C)
{
  static_assert
    (containers::predicates::are_scalar_compatible<
     ::pressio::containers::Tensor<3, A_type>,
     ::pressio::containers::Tensor<3, B_type>,
     ::pressio::containers::Tensor<3, C_type>
     >::value, "Types are not scalar compatible");

  pybind11::object pyblas_ = pybind11::module::import("scipy.linalg.blas");
  const auto & nativeA = *A.data();
  const auto & nativeB = *B.data();
  auto & nativeC = *C.data();

  // if C and B have shape[1] == 1, then we do gemv
  // by slicing A along third axis and
  if (C.extent(1) == 1 and B.extent(1) == 1)
  {
    constexpr auto izero	    = ::pressio::utils::constants<int>::zero();
    constexpr auto ione	    = ::pressio::utils::constants<int>::one();
    constexpr auto transA	    = izero;
    constexpr auto overWritey = ione;
    for (int k=0; k<A.extent(2); ++k)
    {
      // slice A
      const auto tp = pybind11::make_tuple(pybind11::slice(0, A.extent(0), 1),
					   pybind11::slice(0, A.extent(1), 1),
					   k);
      pybind11::array Aslice = nativeA[tp];

      // slice B, C
      const auto tpB = pybind11::make_tuple(pybind11::slice(0, B.extent(0), 1), 0, k);
      pybind11::array BSlice = nativeB[tpB];

      const auto tpC = pybind11::make_tuple(pybind11::slice(0, C.extent(0), 1), 0, k);
      pybind11::array CSlice = nativeC[tpC];

      pyblas_.attr("dgemv")(alpha, Aslice, BSlice, beta, CSlice,
			    izero, ione, izero, ione, transA, overWritey);
    }
  }
  else{
    throw std::runtime_error("Missing impl for extent(1)!=1");
  }
}


/*
 * C = beta * C + alpha*A*B
 * A = rank-3 tensor
 * B,C = rank-2 tensor
*/
template <
  typename A_type, typename B_type, typename scalar_type, typename C_type
  >
::pressio::mpl::enable_if_t<
  containers::predicates::is_fstyle_array_pybind<A_type>::value and
  containers::predicates::is_fstyle_array_pybind<B_type>::value and
  containers::predicates::is_fstyle_array_pybind<C_type>::value
  >
product(::pressio::nontranspose,
	::pressio::nontranspose,
	const scalar_type alpha,
	const ::pressio::containers::Tensor<3, A_type> & A,
	const ::pressio::containers::Tensor<2, B_type> & B,
	const scalar_type beta,
	::pressio::containers::Tensor<2, C_type> & C)
{
  static_assert
    (containers::predicates::are_scalar_compatible<
     ::pressio::containers::Tensor<3, A_type>,
     ::pressio::containers::Tensor<2, B_type>,
     ::pressio::containers::Tensor<2, C_type>
     >::value, "Types are not scalar compatible");

  pybind11::object pyblas_ = pybind11::module::import("scipy.linalg.blas");
  const auto & nativeA = *A.data();
  const auto & nativeB = *B.data();
  auto & nativeC = *C.data();

  constexpr auto izero	    = ::pressio::utils::constants<int>::zero();
  constexpr auto ione	    = ::pressio::utils::constants<int>::one();
  constexpr auto transA	    = izero;
  constexpr auto overWritey = ione;
  for (int k=0; k<A.extent(2); ++k)
  {
    // slice A
    const auto tp = pybind11::make_tuple(pybind11::slice(0, A.extent(0), 1),
					 pybind11::slice(0, A.extent(1), 1),
					 k);
    pybind11::array Aslice = nativeA[tp];

    // slice B
    const auto tpB = pybind11::make_tuple(pybind11::slice(0, B.extent(0), 1), k);
    pybind11::array BSlice = nativeB[tpB];

    // slice C
    const auto tpC = pybind11::make_tuple(pybind11::slice(0, C.extent(0), 1), k);
    pybind11::array CSlice = nativeC[tpC];

    pyblas_.attr("dgemv")(alpha, Aslice, BSlice, beta, CSlice,
			  izero, ione, izero, ione, transA, overWritey);
  }
}

/*
 * C = beta * C + alpha*A^T*B
 * A,B,C = rank-2 tensors
*/
template <typename A_type, typename B_type, typename scalar_type, typename C_type>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<A_type>::value and
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<B_type>::value and
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<C_type>::value
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const scalar_type alpha,
	const A_type & A,
	const B_type & B,
	const scalar_type beta,
	C_type & C)
{
  static_assert
    (containers::predicates::are_scalar_compatible<A_type, B_type, C_type>::value,
     "Types are not scalar compatible");

  // NOTE: need to check if doing this import is expensive,
  // and assess whether we can use blas directly when we know
  // that objects involved are dense with not strange layout.
  pybind11::object pyblas = pybind11::module::import("scipy.linalg.blas");
  constexpr auto one  = ::pressio::utils::constants<scalar_type>::one();
  constexpr auto no   = ::pressio::utils::constants<int>::zero();
  constexpr auto yes  = ::pressio::utils::constants<int>::one();
  constexpr auto transA = yes;
  constexpr auto transB = no;
  constexpr auto ovw    = yes;
  pyblas.attr("dgemm")(one, *A.data(), *B.data(), beta,
		       *C.data(), transA, transB, ovw);
}

/*
 * C = beta * C + alpha*A^T*B
 * A,B = rank-2 tensors
 * C = subspan expression
*/
template <typename A_type, typename B_type, typename scalar_type, typename C_type>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<A_type>::value and
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<B_type>::value
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const scalar_type alpha,
	const A_type & A,
	const B_type & B,
	const scalar_type beta,
	::pressio::containers::expressions::SubspanExpr<C_type> & C)
{
  static_assert
    (containers::predicates::are_scalar_compatible<A_type, B_type, C_type>::value,
     "Types are not scalar compatible");

  for (std::size_t i=0; i<A.extent(1); i++){
    for (std::size_t j=0; j<B.extent(1); j++){
      C(i,j) = beta * C(i,j);
      for (std::size_t k=0; k<A.extent(0); ++k){
	C(i,j) += alpha * A(k,i) * B(k,j);
      }
    }
  }
}


/*
 * C = beta * C + alpha*A^T*B
 * A = rank-3 tensor
 * B,C = rank-2 tensor
*/
template <
  typename A_type, typename B_type, typename scalar_type, typename C_type
  >
::pressio::mpl::enable_if_t<
  containers::predicates::is_fstyle_array_pybind<A_type>::value and
  containers::predicates::is_fstyle_array_pybind<B_type>::value and
  containers::predicates::is_fstyle_array_pybind<C_type>::value
  >
product(::pressio::transpose,
	::pressio::nontranspose,
	const scalar_type alpha,
	const ::pressio::containers::Tensor<3, A_type> & A,
	const ::pressio::containers::Tensor<2, B_type> & B,
	const scalar_type beta,
	::pressio::containers::Tensor<2, C_type> & C)
{
  static_assert
    (containers::predicates::are_scalar_compatible<
     ::pressio::containers::Tensor<3, A_type>,
     ::pressio::containers::Tensor<2, B_type>,
     ::pressio::containers::Tensor<2, C_type>
     >::value, "Types are not scalar compatible");

  pybind11::object pyblas_ = pybind11::module::import("scipy.linalg.blas");

  const auto & nativeA = *A.data();
  const auto & nativeB = *B.data();
  auto & nativeC = *C.data();

  constexpr auto izero	    = ::pressio::utils::constants<int>::zero();
  constexpr auto ione	    = ::pressio::utils::constants<int>::one();
  constexpr auto transA	    = ione;
  constexpr auto overWritey = ione;

  if ( (A.extent(2) == B.extent(1)) and (A.extent(2) == C.extent(1)) )
  {
    for (int k=0; k<A.extent(2); ++k)
      {
	// slice A
	const auto tp = pybind11::make_tuple(pybind11::slice(0, A.extent(0), 1),
					     pybind11::slice(0, A.extent(1), 1),
					     k);
	pybind11::array Aslice = nativeA[tp];

	// slice B
	const auto tpB = pybind11::make_tuple(pybind11::slice(0, B.extent(0), 1), k);
	pybind11::array BSlice = nativeB[tpB];

	// slice C
	const auto tpC = pybind11::make_tuple(pybind11::slice(0, C.extent(0), 1), k);
	pybind11::array CSlice = nativeC[tpC];

	pyblas_.attr("dgemv")(alpha, Aslice, BSlice, beta, CSlice,
			      izero, ione, izero, ione, transA, overWritey);
      }
  }
  else{
    throw std::runtime_error
      ("level-3 product subcase not implemented yet");
  }
}


/***********************************
 * special case A==B and op(A) = transpose
 * C = beta * C + alpha*A^T*A
 **********************************/
template <class A_type, class scalar_type, class C_type>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<A_type>::value and
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<C_type>::value
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const scalar_type alpha,
	const A_type & A,
	const scalar_type beta,
	C_type & C)
{
  // currently not working for expressions because expressions
  // do not have a .data() method and might have non-contiguous layout
  // so we cannot just pass them to blas/lapack, we need to handle them separatly
  static_assert
    (!containers::predicates::is_expression<A_type>::value and
     !containers::predicates::is_expression<C_type>::value,
     "Cannot yet handle expressions for ops::product for pybind11");

  static_assert
    (containers::predicates::are_scalar_compatible<A_type, C_type>::value,
     "Types are not scalar compatible");

  // // attempt to use eigen to map data and do operation but still needs debuggin
  // using mat_t = Eigen::Matrix<scalar_type, -1, -1, Eigen::ColMajor>;
  // Eigen::Map<const mat_t> Am(A.data()->data(), A.extent(0), A.extent(1));
  // auto AmT = Am.transpose();
  // Eigen::Map<mat_t> Cm(C.data()->mutable_data(), C.extent(0), C.extent(1));
  // Cm = beta * Cm + alpha * AmT * Am;

  // NOTE: need to check if doing this import is expensive,
  // and assess whether we can use blas directly when we know
  // that objects involved are dense with not strange layout.
  pybind11::object pyblas = pybind11::module::import("scipy.linalg.blas");
  constexpr auto one  = ::pressio::utils::constants<scalar_type>::one();
  constexpr auto no   = ::pressio::utils::constants<int>::zero();
  constexpr auto yes  = ::pressio::utils::constants<int>::one();
  constexpr auto transA = yes;
  constexpr auto transB = no;
  constexpr auto ovw    = yes;
  pyblas.attr("dgemm")(one, *A.data(), *A.data(), beta,
		       *C.data(), transA, transB, ovw);
}

template <class C_type, class A_type, class scalar_type>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<A_type>::value and
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<C_type>::value,
  C_type
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const scalar_type alpha,
	const A_type & A)
{
  constexpr auto zero  = ::pressio::utils::constants<scalar_type>::zero();
  C_type C(A.extent(1), A.extent(1));
  product(modeA, modeB, alpha, A, zero, C);
  return C;
}

//-------------------------------------------
// C = beta * C + alpha*A*B
// specialize for when A = asDiagonalMatrix expression
//-------------------------------------------
template <class T, class B_type, class scalar_type, class C_type>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_rank1_tensor_wrapper_pybind<T>::value and
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<B_type>::value and
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<C_type>::value
  >
product(::pressio::nontranspose modeA,
	::pressio::nontranspose modeB,
	const scalar_type alpha,
	const pressio::containers::expressions::AsDiagonalMatrixExpr<T> & A,
	const B_type & B,
	const scalar_type beta,
	C_type & C)
{
  static_assert
    (containers::predicates::are_scalar_compatible<T, B_type, C_type>::value,
     "Types are not scalar compatible");

  assert( A.extent(0) == A.extent(1) );
  assert( C.extent(0) == A.extent(0) );
  assert( C.extent(1) == B.extent(1) );
  assert( A.extent(1) == B.extent(0) );
  // since A = asDiagoanlMatrix, get the underlying pressio::Vector
  const auto & Av = *A.pressioObj();

  // we need to make this more efficient, but not now
  for (std::size_t i=0; i<C.extent(0); ++i){
    const auto Avalue = Av(i);
    for (std::size_t j=0; j<C.extent(1); ++j)
    {
      C(i,j) = beta*C(i,j) + alpha*Avalue*B(i,j);
    }
  }
}

}}//end namespace pressio::ops
#endif  // OPS_PYBIND11_OPS_LEVEL3_HPP_
