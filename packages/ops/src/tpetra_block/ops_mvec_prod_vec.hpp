/*
//@HEADER
// ************************************************************************
//
// ops_mvec_prod_vec.hpp
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

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#ifndef OPS_SRC_OPS_TPETRA_BLOCK_MULTI_VECTOR_PROD_VECTOR_HPP_
#define OPS_SRC_OPS_TPETRA_BLOCK_MULTI_VECTOR_PROD_VECTOR_HPP_

namespace pressio{ namespace ops{

/*
 * multi_vector prod vector
 *
 * y = beta * y + alpha*op(A)*x
 *
*/

/* -------------------------------------------------------------------
 * op(A) = A
 * x is a sharedmem vector kokkos wrapper
 *-------------------------------------------------------------------*/
template < typename A_type, typename x_type, typename scalar_type, typename y_type>
::pressio::mpl::enable_if_t<
  containers::meta::is_multi_vector_wrapper_tpetra_block<A_type>::value and
  containers::meta::is_vector_wrapper_kokkos<x_type>::value and
  containers::meta::is_vector_wrapper_tpetra_block<y_type>::value
  >
product(::pressio::nontranspose mode,
	const scalar_type alpha,
	const A_type & A,
	const x_type & x,
	const scalar_type beta,
	y_type & y)
{
  static_assert(containers::meta::are_scalar_compatible<A_type, x_type, y_type>::value,
    "Types are not scalar compatible");

  // make sure the tpetra mv has same exe space of the kokkos vector wrapper
  using tpetra_mvb_dev_t = typename ::pressio::containers::details::traits<A_type>::device_t;
  using kokkos_v_dev_t  = typename ::pressio::containers::details::traits<x_type>::device_type;
  static_assert( std::is_same<tpetra_mvb_dev_t, kokkos_v_dev_t>::value,
  		 "product: tpetra MV and kokkos wrapper need to have same device type" );

  assert( A.numVectors() == x.extent(0) );
  const char ctA = 'N';

  // the the underlying tpetra multivector
  const auto mvViewA = A.data()->getMultiVectorView();

  // get a local view
  const auto yView = y.data()->getVectorView();
  const auto ALocalView_d = mvViewA.getLocalViewDevice();

  // Tpetra::Vector is implemented as a special case of MultiVector //
  // so getLocalView returns a rank-2 view so in order to get
  // view with rank==1 I need to explicitly get the subview of that
  const auto yLocalView_drank2 = yView.getLocalViewDevice();
  const auto yLocalView_drank1 = Kokkos::subview(yLocalView_drank2, Kokkos::ALL(), 0);
  KokkosBlas::gemv(&ctA, alpha, ALocalView_d, *x.data(), beta, yLocalView_drank1);
}


/* -------------------------------------------------------------------
 * op(A) = A
 * x is a sharedmem vector but NOT kokkos
 *-------------------------------------------------------------------*/
template < typename A_type, typename x_type, typename scalar_type, typename y_type>
::pressio::mpl::enable_if_t<
  containers::meta::is_multi_vector_wrapper_tpetra_block<A_type>::value and
  !containers::meta::is_vector_wrapper_kokkos<x_type>::value and
  containers::meta::is_vector_wrapper_tpetra_block<y_type>::value
  >
product(::pressio::nontranspose mode,
	const scalar_type alpha,
	const A_type & A,
	const ::pressio::containers::VectorSharedMemBase<x_type> & x,
	const scalar_type beta,
	y_type & y)
{
  static_assert(containers::meta::are_scalar_compatible<A_type, x_type, y_type>::value,
    "Types are not scalar compatible");

  const auto numVecs = A.numVectors();
  assert(size_t(numVecs) == size_t(x.extent(0)));
  const auto myNrows = A.extentLocal(0);

  // get the wrapped trilinos tpetra multivector
  auto trilD = A.data()->getMultiVectorView();
  auto mv2d = trilD.getLocalViewHost();

  // get wrapped data for the result
  auto y1 = y.data()->getVectorView().getLocalViewHost();
  auto y2 = Kokkos::subview(y1, Kokkos::ALL(), 0);
  y.data()->template modify<Kokkos::HostSpace>();

  // loop
  for (size_t i=0; i<(size_t)myNrows; i++){
    y2[i] = beta*y2[i];
    for (size_t j=0; j<(size_t)numVecs; j++){
      y2[i] += alpha * mv2d(i,j) * x[j];
    }
  }
  using device_t = typename ::pressio::containers::details::traits<y_type>::device_t;
  y.data()->template sync<device_t>();
}


/* -------------------------------------------------------------------
 * x is a distributed Tpetra block vector wrapper
 *-------------------------------------------------------------------*/

// y = sharedmem vec not kokkos
template <typename A_type, typename x_type, typename y_type, typename scalar_type>
::pressio::mpl::enable_if_t<
  containers::meta::is_multi_vector_wrapper_tpetra_block<A_type>::value and
  containers::meta::is_vector_wrapper_tpetra_block<x_type>::value and
  !containers::meta::is_vector_wrapper_kokkos<y_type>::value
  >
product(::pressio::transpose mode,
	const scalar_type alpha,
	const A_type & A,
	const x_type & x,
	const scalar_type beta,
	::pressio::containers::VectorSharedMemBase<y_type> & y)
{
  /* workaround the non-constness of getVectorView*/
  using wrapped_t = typename containers::details::traits<x_type>::wrapped_t;
  const auto xvv = const_cast<wrapped_t*>(x.data())->getVectorView();
  const auto mvA_mvv = A.data()->getMultiVectorView();
  const auto numVecs = A.numVectors();
  for (std::size_t i=0; i<(std::size_t)numVecs; i++){
    // colI is a Teuchos::RCP<Vector<...>>
    const auto colI = mvA_mvv.getVector(i);
    y[i] = beta*y[i] + alpha * colI->dot(xvv);
  }
}

// y = wrapper of Kokkos vector
template <typename A_type, typename x_type, typename y_type, typename scalar_type>
::pressio::mpl::enable_if_t<
  containers::meta::is_multi_vector_wrapper_tpetra_block<A_type>::value and
  containers::meta::is_vector_wrapper_tpetra_block<x_type>::value and
  containers::meta::is_vector_wrapper_kokkos<y_type>::value
  >
product(::pressio::transpose mode,
	const scalar_type alpha,
	const A_type & A,
	const x_type & x,
	const scalar_type beta,
	y_type & y)
{
  constexpr auto zero = ::pressio::utils::constants::zero<scalar_type>();
  constexpr auto one  = ::pressio::utils::constants::one<scalar_type>();
  if (alpha != one and beta != zero)
    throw std::runtime_error("y = beta * y + alpha*op(A)*x where A = Tpetra MV wrapper and \
x=Tpetra Vector wrapper and y=Kokkos wrapper is not yet supported for beta!=0 and alpha!=1.");

  const auto A_mvv = A.data()->getMultiVectorView();
  using tpetra_blockvector_t = typename containers::details::traits<x_type>::wrapped_t;
  const auto x_vv = const_cast<tpetra_blockvector_t*>(x.data())->getVectorView();
  auto request = Tpetra::idot( *y.data(), A_mvv, x_vv);
  request->wait();
}

}}//end namespace pressio::ops
#endif
#endif
