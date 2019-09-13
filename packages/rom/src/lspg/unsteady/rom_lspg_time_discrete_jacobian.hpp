/*
//@HEADER
// ************************************************************************
//
// rom_lspg_time_discrete_jacobian.hpp
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

#ifndef ROM_LSPG_TIME_DISCRETE_JACOBIAN_HPP_
#define ROM_LSPG_TIME_DISCRETE_JACOBIAN_HPP_

#include "../../../../ode/src/ode_ConfigDefs.hpp"
#include "../../../../ode/src/implicit/ode_implicit_constants.hpp"

namespace pressio{ namespace rom{ namespace impl{

template <
  ode::ImplicitEnum odeMethod,
  typename lspg_matrix_type,
  typename scalar_type,
  typename decoder_jac_type,
  typename ud_ops
#ifdef HAVE_PYBIND11
  , mpl::enable_if_t<
    mpl::not_same< ud_ops, pybind11::object>::value
    > * = nullptr
#endif
  >
void time_discrete_jacobian(lspg_matrix_type & jphi, //jphi holds J * phi
			    scalar_type	dt,
			    const decoder_jac_type & phi,
			    const ud_ops * udOps){

  // prefactor (f) multiplying f*dt*J*phi
  auto prefactor = static_cast<scalar_type>(1);
  if (odeMethod == ode::ImplicitEnum::BDF2)
    prefactor = ode::coeffs::bdf2<scalar_type>::c3_;

  udOps->time_discrete_jacobian(*jphi.data(), *phi.data(), prefactor, dt);
}


#ifdef HAVE_PYBIND11
template <
  ode::ImplicitEnum odeMethod,
  typename lspg_matrix_type,
  typename scalar_type,
  typename decoder_jac_type,
  typename ud_ops,
  mpl::enable_if_t<
    mpl::is_same< ud_ops, pybind11::object>::value
    > * = nullptr
  >
void time_discrete_jacobian(lspg_matrix_type & jphi, //jphi holds J * phi
			    scalar_type	dt,
			    const decoder_jac_type & phi,
			    const ud_ops & udOps){

  // prefactor (f) multiplying f*dt*J*phi
  auto prefactor = static_cast<scalar_type>(1);
  if (odeMethod == ode::ImplicitEnum::BDF2)
    prefactor = ode::coeffs::bdf2<scalar_type>::c3_;

  udOps.attr("time_discrete_jacobian")(jphi, phi, prefactor, dt);
}
#endif



/*
 * for EIGEN
 * only if jphi, phi are of same size
*/
template <
  ode::ImplicitEnum odeMethod,
  typename lspg_matrix_type,
  typename scalar_type,
  typename decoder_jac_type,
  ::pressio::mpl::enable_if_t<
    (containers::meta::is_multi_vector_wrapper_eigen<lspg_matrix_type>::value and
     containers::meta::is_multi_vector_wrapper_eigen<decoder_jac_type>::value)
#ifdef HAVE_KOKKOS
    or
    (containers::meta::is_multi_vector_wrapper_kokkos<lspg_matrix_type>::value and
     containers::meta::is_multi_vector_wrapper_kokkos<decoder_jac_type>::value)
#endif
    > * = nullptr
  >
void time_discrete_jacobian(lspg_matrix_type & jphi, //jphi holds J * phi
			    scalar_type	dt,
			    const decoder_jac_type & phi){

  assert( jphi.numVectors() == phi.numVectors() );
  assert( jphi.length() == phi.length() );

  // prefactor (f) multiplying f*dt*J*phi
  auto prefactor = ::pressio::utils::constants::one<scalar_type>();
  if (odeMethod == ode::ImplicitEnum::BDF2)
    prefactor = ode::coeffs::bdf2<scalar_type>::c3_;

  constexpr auto one = ::pressio::utils::constants::one<scalar_type>();
  constexpr auto negone = ::pressio::utils::constants::negOne<scalar_type>();
  const auto a = negone*prefactor*dt;

  //loop over elements of jphi
  ::pressio::containers::ops::do_update(jphi, a, phi, one);
  // for (auto i=0; i<jphi.length(); i++){
  //   for (auto j=0; j<jphi.numVectors(); j++)
  //     jphi(i,j) = phi(i,j) - prefactor*dt*jphi(i,j);
  // }
}



/* when we have hyper-reduction and we need to calculate the
 * time-discrete residual, yn (i.e. the current state) and
 * ynm (i.e. the states at prev steps) might have
 * different entries/maps than R (i.e. the spatial residual)
 * so we need to update R only picking elements corresponding
 * elements in yn and ynm.
 * In other words, in the functions below we do NOT assume that
 * yn and ynm have the same distribution/sizes of R but we DO assume
 * that yn and ynm contain at least all the elements of R (and more).
 * THis is because to evaluate spatial residual one typically
 * needs neighbors values, so in the hyper-red scenario the
 * states contain more DOFs than the spatial residual.
 * for instance, imagine the i-th cell in a uniform grid with a
 * centered 2nd-order FD operator.
 * To compute the spatial residual at the i-th cell, one needs the
 * cell value y_i, as well as the left, y_i-1, and right, y_i+1,
 * values of the state. So for a single DOF of the residual, one
 * needs 3 DOFs of the state. Extrapolating to all cells, then we
 * clearly see that we need more dofs for the states than for the residual.
 *
 * Below, we account for this potential scenario and the code below
 * obviously also works for the case where the vectors have same
 * the same distributions.
*/

#ifdef HAVE_TRILINOS


/*************************************
            epetra
*************************************/
template <
  ode::ImplicitEnum odeMethod,
  typename lspg_matrix_type,
  typename scalar_type,
  typename decoder_jac_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_epetra<lspg_matrix_type>::value and
    containers::meta::is_multi_vector_wrapper_epetra<decoder_jac_type>::value
    > * = nullptr
  >
void time_discrete_jacobian(lspg_matrix_type & jphi, //jphi stands for J * phi
			    scalar_type	dt,
			    const decoder_jac_type & phi){

  // integral type of the global indices
  using GO_t = typename containers::details::traits<lspg_matrix_type>::global_ordinal_t;

  //get row map of phi
  const auto & phi_map = phi.getDataMap();
  // get my global elements
  std::vector<GO_t> gIDphi( phi.localLength() );
  phi_map.MyGlobalElements( gIDphi.data() );

  // get map of jphi
  const auto & jphi_map = jphi.getDataMap();
  // get global elements
  std::vector<GO_t> gIDjphi( jphi.localLength() );
  jphi_map.MyGlobalElements( gIDjphi.data() );

  // prefactor (f) multiplying f*dt*J*phi
  auto prefactor = static_cast<scalar_type>(1);
  if (odeMethod == ode::ImplicitEnum::BDF2)
    prefactor = ode::coeffs::bdf2<scalar_type>::c3_;

  //loop over elements of jphi
  for (auto i=0; i<jphi.localLength(); i++){
    // ask the phi map what is the local index corresponding
    // to the global index we are handling
    auto lid = phi_map.LID(gIDjphi[i]);
    for (auto j=0; j<jphi.globalNumVectors(); j++)
      jphi(i,j) = phi(lid,j) - prefactor*dt*jphi(i,j);
  }
}

/*************************************
            tpetra
*************************************/
template <
  ode::ImplicitEnum odeMethod,
  typename lspg_matrix_type,
  typename scalar_type,
  typename decoder_jac_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_tpetra<lspg_matrix_type>::value and
    containers::meta::is_multi_vector_tpetra<decoder_jac_type>::value
    > * = nullptr
  >
void time_discrete_jacobian(lspg_matrix_type & jphi, //jphi holds J * phi
			    scalar_type	dt,
			    const decoder_jac_type & phi){

  //get row map of phi
  const auto phi_map = phi.getMap();
  // get my global elements
  const auto gIDphi = phi_map->getMyGlobalIndices();

  // get map of jphi
  const auto jphi_map = jphi.getMap();
  // get global elements
  const auto gIDjphi = jphi_map->getMyGlobalIndices();

  // prefactor (f) multiplying f*dt*J*phi
  auto prefactor = static_cast<scalar_type>(1);
  if (odeMethod == ode::ImplicitEnum::BDF2)
    prefactor = ode::coeffs::bdf2<scalar_type>::c3_;

  auto jphi2dView = jphi.get2dViewNonConst();
  auto phi2dView = phi.get2dView();

  //loop over elements of jphi
  for (size_t i=0; i<jphi.getLocalLength(); i++){
    // ask the phi map what is the local index corresponding
    // to the global index we are handling
    const auto lid = phi_map->getLocalElement(gIDjphi[i]);
    for (size_t j=0; j<jphi.getNumVectors(); j++)
      jphi2dView[j][i] = phi2dView[j][lid] - prefactor*dt*jphi2dView[j][i];
  }
}

/*************************************
            tpetra block
*************************************/
template <
  ode::ImplicitEnum odeMethod,
  typename lspg_matrix_type,
  typename scalar_type,
  typename decoder_jac_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_tpetra<lspg_matrix_type>::value and
    containers::meta::is_multi_vector_wrapper_tpetra<decoder_jac_type>::value
    > * = nullptr
  >
void time_discrete_jacobian(lspg_matrix_type & jphi,
			    scalar_type	dt,
			    const decoder_jac_type & phi){

  time_discrete_jacobian<odeMethod>(*jphi.data(), dt, *phi.data());
}

template <
  ode::ImplicitEnum odeMethod,
  typename lspg_matrix_type,
  typename scalar_type,
  typename decoder_jac_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_tpetra_block<lspg_matrix_type>::value and
    containers::meta::is_multi_vector_wrapper_tpetra_block<decoder_jac_type>::value
    > * = nullptr
  >
void time_discrete_jacobian(lspg_matrix_type & jphi,
			    scalar_type	dt,
			    const decoder_jac_type & phi){
  auto jphi_mvv = jphi.data()->getMultiVectorView();
  auto phi_mvv  = phi.data()->getMultiVectorView();
  time_discrete_jacobian<odeMethod>(jphi_mvv, dt, phi_mvv);
}

#endif

}}}//end namespace pressio::rom::impl
#endif
