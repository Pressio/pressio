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

#ifndef ROM_LSPG_IMPL_UNSTEADY_CONTINUOUS_TIME_API_DISCRETE_TIME_FUNCTIONS_ROM_LSPG_TIME_DISCRETE_JACOBIAN_HPP_
#define ROM_LSPG_IMPL_UNSTEADY_CONTINUOUS_TIME_API_DISCRETE_TIME_FUNCTIONS_ROM_LSPG_TIME_DISCRETE_JACOBIAN_HPP_

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{ namespace unsteady{

template <typename stepper_tag, typename scalar_t> struct dtPrefactor;

template <typename scalar_t>
struct dtPrefactor<::pressio::ode::implicitmethods::Euler, scalar_t>{
  static constexpr auto value = ::pressio::ode::constants::bdf1<scalar_t>::c_f_;
};

template <typename scalar_t>
struct dtPrefactor<::pressio::ode::implicitmethods::BDF2, scalar_t>{
  static constexpr auto value = ::pressio::ode::constants::bdf2<scalar_t>::c_f_;
};

template <typename scalar_t>
struct dtPrefactor<::pressio::ode::implicitmethods::CrankNicolson, scalar_t>{
  static constexpr auto value = ::pressio::ode::constants::cranknicolson<scalar_t>::c_fnp1_;
};
// ------------------------------------------------------

// regular c++ with user-defined OPS
template <
  typename stepper_tag,
  typename lspg_matrix_type,
  typename scalar_type,
  typename decoder_jac_type,
  typename ud_ops
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  , mpl::enable_if_t<
     !::pressio::containers::predicates::is_tensor_wrapper_pybind<lspg_matrix_type>::value and
     mpl::not_same< ud_ops, pybind11::object>::value, int > = 0
#endif
  >
void time_discrete_jacobian(lspg_matrix_type & jphi, //jphi holds J * phi
			    const scalar_type	& dt,
			    const decoder_jac_type & phi,
			    const ud_ops * udOps)
{

  constexpr auto one = ::pressio::utils::constants<scalar_type>::one();
  const auto prefactor = dt * dtPrefactor<stepper_tag, scalar_type>::value;
  udOps->time_discrete_jacobian(prefactor, *jphi.data(), one, *phi.data());
}


#ifdef PRESSIO_ENABLE_TPL_PYBIND11
template <
  typename stepper_tag,
  typename lspg_matrix_type,
  typename scalar_type,
  typename decoder_jac_type
  >
mpl::enable_if_t<
 ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<lspg_matrix_type>::value
  >
time_discrete_jacobian(lspg_matrix_type & jphi, //jphi holds J * phi
		       const scalar_type	& dt,
		       const decoder_jac_type & phi)
{
  // prefactor (f) multiplying f*dt*J*phi
  const auto prefactor = dt * dtPrefactor<stepper_tag, scalar_type>::value;
  const auto nRows = jphi.extent(0);
  const auto nCols = jphi.extent(1);
  for (std::size_t j=0; j<(std::size_t)nCols; ++j){
    for (std::size_t i=0; i<(std::size_t)nRows; ++i){
      jphi(i,j) = phi(i,j) + prefactor*jphi(i,j);
    }
  }
}

template <
  typename stepper_tag,
  typename lspg_matrix_type,
  typename scalar_type,
  typename decoder_jac_type,
  typename hyp_ind_t
  >
mpl::enable_if_t<
  ::pressio::containers::predicates::is_fstyle_rank2_tensor_wrapper_pybind<lspg_matrix_type>::value
  >
time_discrete_jacobian(lspg_matrix_type & jphi, //jphi holds J * phi
		       const scalar_type	& dt,
		       const decoder_jac_type & phi,
		       const hyp_ind_t & hypIndices)
{
  // hypindices has same extent as sample mesh and contains
  // indices of the entries in the state that correspond to
  // the sample mesh points
  assert(jphi.extent(0) == hypIndices.extent(0));

  const auto prefactor = dt * dtPrefactor<stepper_tag, scalar_type>::value;
  const auto nRows = jphi.extent(0);
  const auto nCols = jphi.extent(1);
  for (std::size_t j=0; j<(std::size_t)nCols; ++j)
  {
    for (std::size_t i=0; i<(std::size_t)nRows; ++i)
    {
      const auto rowInd = hypIndices(i);
      jphi(i,j) = phi(rowInd,j) + prefactor*jphi(i,j);
    }
  }
}
#endif


#ifdef PRESSIO_ENABLE_TPL_EIGEN
template <
  typename stepper_tag,
  typename lspg_matrix_type,
  typename scalar_type,
  typename decoder_jac_type,
  typename hyp_ind_t
  >
mpl::enable_if_t<
  ::pressio::containers::predicates::is_dense_matrix_wrapper_eigen<lspg_matrix_type>::value or
  ::pressio::containers::predicates::is_multi_vector_wrapper_eigen<lspg_matrix_type>::value
  >
time_discrete_jacobian(lspg_matrix_type & jphi, //jphi holds J * phi
		       const scalar_type	& dt,
		       const decoder_jac_type & phi,
		       const hyp_ind_t & hypIndices)
{
  // hypindices has same extent as sample mesh and contains
  // indices of the entries in the state that correspond to
  // the sample mesh points
  assert(jphi.extent(0) == hypIndices.extent(0));

  const auto prefactor = dt * dtPrefactor<stepper_tag, scalar_type>::value;
  const auto nRows = jphi.extent(0);
  const auto nCols = jphi.extent(1);
  for (std::size_t i=0; i<(std::size_t)nRows; ++i)
  {
    const auto rowInd = hypIndices(i);
    for (std::size_t j=0; j<(std::size_t)nCols; ++j){
      jphi(i,j) = phi(rowInd,j) + prefactor*jphi(i,j);
    }
  }
}
#endif

/*
 * for EIGEN or Kokkos, with Euler and no hyper-reduction
 * only if jphi, phi are of same size
*/
#ifdef PRESSIO_ENABLE_TPL_EIGEN
template <
  typename stepper_tag,
  typename lspg_matrix_type,
  typename scalar_type,
  typename decoder_jac_type
>
::pressio::mpl::enable_if_t<
  containers::predicates::is_multi_vector_wrapper_eigen<lspg_matrix_type>::value and
  containers::predicates::is_multi_vector_wrapper_eigen<decoder_jac_type>::value
>
time_discrete_jacobian(lspg_matrix_type & jphi, //jphi holds J * phi
		       const scalar_type	& dt,
		       const decoder_jac_type & phi)
{

  assert( jphi.numVectors() == phi.numVectors() );
  assert( jphi.extent(0) == phi.extent(0) );

  // prefactor (f) multiplying f*dt*J*phi
  const auto prefactor = dt * dtPrefactor<stepper_tag, scalar_type>::value;

  // jphi = phi + prefactor*dt*jphi
  constexpr auto one = ::pressio::utils::constants<scalar_type>::one();
  ::pressio::ops::update(jphi, prefactor, phi, one);
}
#endif

#ifdef PRESSIO_ENABLE_TPL_KOKKOS
template <
  typename stepper_tag,
  typename lspg_matrix_type,
  typename scalar_type,
  typename decoder_jac_type
>
::pressio::mpl::enable_if_t<
  containers::predicates::is_multi_vector_wrapper_kokkos<lspg_matrix_type>::value and
  containers::predicates::is_multi_vector_wrapper_kokkos<decoder_jac_type>::value
>
time_discrete_jacobian(lspg_matrix_type & jphi, //jphi holds J * phi
		       const scalar_type	& dt,
		       const decoder_jac_type & phi)
{

  assert( jphi.numVectors() == phi.numVectors() );
  assert( jphi.extent(0) == phi.extent(0) );

  // prefactor (f) multiplying f*dt*J*phi
  const auto prefactor = dt * dtPrefactor<stepper_tag, scalar_type>::value;

  // jphi = phi + prefactor*dt*jphi
  constexpr auto one = ::pressio::utils::constants<scalar_type>::one();
  ::pressio::ops::update(jphi, prefactor, phi, one);
}
#endif


#ifdef PRESSIO_ENABLE_TPL_TRILINOS

/* for trilinos data structures, support hyper-reduction too.
 * When we have hyper-reduction and we need to calculate the
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


/*************************************
            epetra
*************************************/
template <
  typename stepper_tag,
  typename lspg_matrix_type,
  typename scalar_type,
  typename decoder_jac_type
>
::pressio::mpl::enable_if_t<
  containers::predicates::is_multi_vector_wrapper_epetra<lspg_matrix_type>::value and
  containers::predicates::is_multi_vector_wrapper_epetra<decoder_jac_type>::value
  >
time_discrete_jacobian(lspg_matrix_type & jphi, //jphi stands for J * phi
		       const scalar_type & dt,
		       const decoder_jac_type & phi)
{

  // integral type of the global indices
  using GO_t = typename lspg_matrix_type::traits::global_ordinal_t;

  // row map of phi
  const auto & phi_map = phi.data()->Map();
  // my global elements
  std::vector<GO_t> gIDphi( phi.extentLocal(0) );
  phi_map.MyGlobalElements( gIDphi.data() );

  // map of jphi
  const auto & jphi_map = jphi.data()->Map();
  // global elements
  std::vector<GO_t> gIDjphi( jphi.extentLocal(0) );
  jphi_map.MyGlobalElements( gIDjphi.data() );

  // prefactor (f) multiplying f*dt*J*phi
  constexpr auto prefactor = dtPrefactor<stepper_tag, scalar_type>::value;

  //loop over elements of jphi
  for (auto i=0; i<jphi.extentLocal(0); i++){
    // ask the phi map what is the local index corresponding
    // to the global index we are handling
    auto lid = phi_map.LID(gIDjphi[i]);
    for (auto j=0; j<jphi.numVectors(); j++)
      jphi(i,j) = phi(lid,j) + prefactor*dt*jphi(i,j);
  }
}

/*************************************
            tpetra
*************************************/
template <
  typename stepper_tag,
  typename lspg_matrix_type,
  typename scalar_type,
  typename decoder_jac_type
>
::pressio::mpl::enable_if_t<
  containers::predicates::is_multi_vector_tpetra<lspg_matrix_type>::value and
  containers::predicates::is_multi_vector_tpetra<decoder_jac_type>::value
  >
time_discrete_jacobian(lspg_matrix_type & jphi, //jphi holds J * phi
		       const scalar_type & dt,
		       const decoder_jac_type & phi)
{

  // row map of phi
  const auto phi_map = phi.getMap();
  // my global elements
  const auto gIDphi = phi_map->getMyGlobalIndices();

  // map of jphi
  const auto jphi_map = jphi.getMap();
  // global elements
  const auto gIDjphi = jphi_map->getMyGlobalIndices();

  // prefactor (f) multiplying f*dt*J*phi
  constexpr auto prefactor = dtPrefactor<stepper_tag, scalar_type>::value;

  auto jphi2dView = jphi.get2dViewNonConst();
  auto phi2dView = phi.get2dView();

  //loop over elements of jphi
  for (size_t i=0; i<jphi.getLocalLength(); i++){
    // ask the phi map what is the local index corresponding
    // to the global index we are handling
    const auto lid = phi_map->getLocalElement(gIDjphi[i]);
    for (size_t j=0; j<jphi.getNumVectors(); j++)
      jphi2dView[j][i] = phi2dView[j][lid] + prefactor*dt*jphi2dView[j][i];
  }
}


/*************************************
            tpetra block
*************************************/
template <
  typename stepper_tag,
  typename lspg_matrix_type,
  typename scalar_type,
  typename decoder_jac_type
>
::pressio::mpl::enable_if_t<
  containers::predicates::is_multi_vector_wrapper_tpetra<lspg_matrix_type>::value and
  containers::predicates::is_multi_vector_wrapper_tpetra<decoder_jac_type>::value
>
time_discrete_jacobian(lspg_matrix_type & jphi,
			    const scalar_type	& dt,
			    const decoder_jac_type & phi)
{
  time_discrete_jacobian<stepper_tag>(*jphi.data(), dt, *phi.data());
}

template <
  typename stepper_tag,
  typename lspg_matrix_type,
  typename scalar_type,
  typename decoder_jac_type
>
::pressio::mpl::enable_if_t<
  containers::predicates::is_multi_vector_wrapper_tpetra_block<lspg_matrix_type>::value and
  containers::predicates::is_multi_vector_wrapper_tpetra_block<decoder_jac_type>::value
>
time_discrete_jacobian(lspg_matrix_type & jphi,
			    const scalar_type & dt,
			    const decoder_jac_type & phi)
{
  auto jphi_mvv = jphi.data()->getMultiVectorView();
  auto phi_mvv  = phi.data()->getMultiVectorView();
  time_discrete_jacobian<stepper_tag>(jphi_mvv, dt, phi_mvv);
}

#endif

}}}}}//end namespace pressio::rom::lspg::unsteady::impl
#endif  // ROM_LSPG_IMPL_UNSTEADY_CONTINUOUS_TIME_API_DISCRETE_TIME_FUNCTIONS_ROM_LSPG_TIME_DISCRETE_JACOBIAN_HPP_
