/*
//@HEADER
// ************************************************************************
//
// rom_lspg_time_discrete_residual.hpp
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

#ifndef ROM_LSPG_IMPL_UNSTEADY_CONTINUOUS_TIME_API_DISCRETE_TIME_FUNCTIONS_ROM_LSPG_TIME_DISCRETE_RESIDUAL_HPP_
#define ROM_LSPG_IMPL_UNSTEADY_CONTINUOUS_TIME_API_DISCRETE_TIME_FUNCTIONS_ROM_LSPG_TIME_DISCRETE_RESIDUAL_HPP_

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{ namespace unsteady{

/* enable when BDF1 and user-defined ops */
template<
  typename stepper_tag,
  typename fom_states_manager_t,
  typename state_type,
  typename scalar_type,
  typename ud_ops
  >
::pressio::mpl::enable_if_t<
  std::is_same<stepper_tag, ::pressio::ode::implicitmethods::Euler>::value and
  mpl::not_void<ud_ops>::value
  >
time_discrete_residual(const fom_states_manager_t & fomStatesMngr,
		       state_type & R,
		       const scalar_type & dt,
		       const ud_ops * udOps)
{

  const auto & y_np1  = fomStatesMngr.fomStateAt(::pressio::ode::nPlusOne());
  const auto & y_n    = fomStatesMngr.fomStateAt(::pressio::ode::n());
  constexpr auto cnp1 = ::pressio::ode::constants::bdf1<scalar_type>::c_np1_;
  constexpr auto cn   = ::pressio::ode::constants::bdf1<scalar_type>::c_n_;
  const auto cfdt     = ::pressio::ode::constants::bdf1<scalar_type>::c_f_ * dt;

  udOps->time_discrete_residual(cfdt, *R.data(),
				cnp1, *y_np1.data(),
				cn, *y_n.data());
}

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
/* default BDF1 */
template<
  typename stepper_tag,
  typename fom_states_manager_t,
  typename residual_type,
  typename scalar_type
  >
::pressio::mpl::enable_if_t<
  std::is_same<stepper_tag, ::pressio::ode::implicitmethods::Euler>::value and
  ::pressio::containers::predicates::is_rank1_tensor_wrapper_pybind<residual_type>::value
  >
time_discrete_residual(const fom_states_manager_t & fomStatesMngr,
		       residual_type & R,
		       const scalar_type & dt)
{
  const auto & y_np1  = fomStatesMngr.fomStateAt(::pressio::ode::nPlusOne());
  const auto & y_n    = fomStatesMngr.fomStateAt(::pressio::ode::n());

  constexpr auto cnp1 = ::pressio::ode::constants::bdf1<scalar_type>::c_np1_;
  constexpr auto cn   = ::pressio::ode::constants::bdf1<scalar_type>::c_n_;
  const auto cfdt     = ::pressio::ode::constants::bdf1<scalar_type>::c_f_ * dt;

  ::pressio::ops::update(R, cfdt, y_np1, cnp1, y_n, cn);
}

/* default BDF2 */
template<
  typename stepper_tag,
  typename fom_states_manager_t,
  typename residual_type,
  typename scalar_type
  >
::pressio::mpl::enable_if_t<
  std::is_same<stepper_tag, ::pressio::ode::implicitmethods::BDF2>::value and
  ::pressio::containers::predicates::is_rank1_tensor_wrapper_pybind<residual_type>::value
  >
time_discrete_residual(const fom_states_manager_t & fomStatesMngr,
		       residual_type & R,
		       const scalar_type & dt)
{
  const auto & y_np1 = fomStatesMngr.fomStateAt(::pressio::ode::nPlusOne());
  const auto & y_n   = fomStatesMngr.fomStateAt(::pressio::ode::n());
  const auto & y_nm1 = fomStatesMngr.fomStateAt(::pressio::ode::nMinusOne());

  constexpr auto cnp1 = ::pressio::ode::constants::bdf2<scalar_type>::c_np1_;
  constexpr auto cn   = ::pressio::ode::constants::bdf2<scalar_type>::c_n_;
  constexpr auto cnm1 = ::pressio::ode::constants::bdf2<scalar_type>::c_nm1_;
  const auto cfdt     = ::pressio::ode::constants::bdf2<scalar_type>::c_f_ * dt;

  ::pressio::ops::update(R, cfdt,
			 y_np1, cnp1,
			 y_n, cn,
			 y_nm1, cnm1);
}

/* BDF1 and hyper-reduction indices passed in */
template<
  typename stepper_tag,
  typename fom_states_manager_t,
  typename residual_type,
  typename scalar_type,
  typename hyp_ind_t
  >
::pressio::mpl::enable_if_t<
  std::is_same<stepper_tag, ::pressio::ode::implicitmethods::Euler>::value and
   ::pressio::containers::predicates::is_rank1_tensor_wrapper_pybind<residual_type>::value
  >
time_discrete_residual(const fom_states_manager_t & fomStatesMngr,
		       residual_type & R,
		       const scalar_type & dt,
		       const hyp_ind_t & hypIndices)
{
  const auto & y_np1 = fomStatesMngr.fomStateAt(::pressio::ode::nPlusOne());
  const auto & y_n   = fomStatesMngr.fomStateAt(::pressio::ode::n());
  constexpr auto cnp1 = ::pressio::ode::constants::bdf1<scalar_type>::c_np1_;
  constexpr auto cn   = ::pressio::ode::constants::bdf1<scalar_type>::c_n_;
  const auto cfdt     = ::pressio::ode::constants::bdf1<scalar_type>::c_f_ * dt;

  // hypindices has same extent as sample mesh and contains
  // indices of the entries in the state that correspond to
  // the sample mesh points

  //R = y_n+1 - y_n - dt * R but we need to combine the correct entries
  assert(R.extent(0) == hypIndices.extent(0));
  for (std::size_t i=0; i<(std::size_t) R.extent(0); ++i)
  {
    const auto yI = hypIndices(i);
    R(i) = cnp1*y_np1(yI) + cn*y_n(yI) + cfdt*R(i);
  }
}

/* BDF2 and hyper-reduction indices passed in */
template<
  typename stepper_tag,
  typename fom_states_manager_t,
  typename residual_type,
  typename scalar_type,
  typename hyp_ind_t
  >
::pressio::mpl::enable_if_t<
  std::is_same<stepper_tag, ::pressio::ode::implicitmethods::BDF2>::value and
  ::pressio::containers::predicates::is_rank1_tensor_wrapper_pybind<residual_type>::value
  >
time_discrete_residual(const fom_states_manager_t & fomStatesMngr,
		       residual_type & R,
		       const scalar_type & dt,
		       const hyp_ind_t & hypIndices)
{
  const auto & y_np1 = fomStatesMngr.fomStateAt(::pressio::ode::nPlusOne());
  const auto & y_n   = fomStatesMngr.fomStateAt(::pressio::ode::n());
  const auto & y_nm1 = fomStatesMngr.fomStateAt(::pressio::ode::nMinusOne());

  constexpr auto cnp1 = ::pressio::ode::constants::bdf2<scalar_type>::c_np1_;
  constexpr auto cn   = ::pressio::ode::constants::bdf2<scalar_type>::c_n_;
  constexpr auto cnm1 = ::pressio::ode::constants::bdf2<scalar_type>::c_nm1_;
  const auto cfdt     = ::pressio::ode::constants::bdf2<scalar_type>::c_f_ * dt;

  // hypindices has same extent as sample mesh and contains
  // indices of the entries in the state that correspond to
  // the sample mesh points
  assert(R.extent(0) == hypIndices.extent(0));
  for (std::size_t i=0; i<(std::size_t) R.extent(0); ++i)
  {
    const auto yI = hypIndices(i);
    R(i) = cnp1*y_np1(yI) + cn*y_n(yI) + cnm1*y_nm1(yI) + cfdt*R(i);
  }
}
#endif


// ----------------------------------------------------------------------
// EIGEN with hyp red enabled with vec of indices
// ----------------------------------------------------------------------
#ifdef PRESSIO_ENABLE_TPL_EIGEN
/* BDF1 and hyper-reduction indices passed in */
template<
  typename stepper_tag,
  typename fom_states_manager_t,
  typename residual_type,
  typename scalar_type,
  typename hyp_ind_t
  >
::pressio::mpl::enable_if_t<
  std::is_same<stepper_tag, ::pressio::ode::implicitmethods::Euler>::value and
  ::pressio::containers::predicates::is_vector_wrapper_eigen<residual_type>::value
  >
time_discrete_residual(const fom_states_manager_t & fomStatesMngr,
		       residual_type & R,
		       const scalar_type & dt,
		       const hyp_ind_t & hypIndices)
{
  const auto & y_np1 = fomStatesMngr.fomStateAt(::pressio::ode::nPlusOne());
  const auto & y_n   = fomStatesMngr.fomStateAt(::pressio::ode::n());

  constexpr auto cnp1 = ::pressio::ode::constants::bdf1<scalar_type>::c_np1_;
  constexpr auto cn   = ::pressio::ode::constants::bdf1<scalar_type>::c_n_;
  const auto cfdt     = ::pressio::ode::constants::bdf1<scalar_type>::c_f_ * dt;

  // hypindices has same extent as sample mesh and contains
  // indices of the entries in the state that correspond to
  // the sample mesh points
  assert(R.extent(0) == hypIndices.extent(0));
  for (std::size_t i=0; i<(std::size_t) R.extent(0); ++i)
  {
    const auto yI = hypIndices(i);
    R(i) = cnp1*y_np1(yI) + cn*y_n(yI) + cfdt*R(i);
  }
}


/* BDF2 and hyper-reduction indices passed in */
template<
  typename stepper_tag,
  typename fom_states_manager_t,
  typename residual_type,
  typename scalar_type,
  typename hyp_ind_t
  >
::pressio::mpl::enable_if_t<
  std::is_same<stepper_tag, ::pressio::ode::implicitmethods::BDF2>::value and
  ::pressio::containers::predicates::is_vector_wrapper_eigen<residual_type>::value
  >
time_discrete_residual(const fom_states_manager_t & fomStatesMngr,
		       residual_type & R,
		       const scalar_type & dt,
		       const hyp_ind_t & hypIndices)
{
  const auto & y_np1 = fomStatesMngr.fomStateAt(::pressio::ode::nPlusOne());
  const auto & y_n   = fomStatesMngr.fomStateAt(::pressio::ode::n());
  const auto & y_nm1 = fomStatesMngr.fomStateAt(::pressio::ode::nMinusOne());

  constexpr auto cnp1 = ::pressio::ode::constants::bdf2<scalar_type>::c_np1_;
  constexpr auto cn   = ::pressio::ode::constants::bdf2<scalar_type>::c_n_;
  constexpr auto cnm1 = ::pressio::ode::constants::bdf2<scalar_type>::c_nm1_;

  const auto cfdt     = ::pressio::ode::constants::bdf2<scalar_type>::c_f_ * dt;

  // hypindices has same extent as sample mesh and contains
  // indices of the entries in the state that correspond to
  // the sample mesh points
  assert(R.extent(0) == hypIndices.extent(0));
  for (std::size_t i=0; i<(std::size_t) R.extent(0); ++i)
  {
    const auto yI = hypIndices(i);
    R(i) = cnp1*y_np1(yI) + cn*y_n(yI) + cnm1*y_nm1(yI) + cfdt*R(i);
  }
}

// Crank-Nicolson and hyper-reduction indices passed in */
template<
  typename stepper_tag,
  typename fom_states_manager_t,
  typename fom_velocities_manager_t,
  typename residual_type,
  typename scalar_type,
  typename hyp_ind_t
  >
::pressio::mpl::enable_if_t<
  std::is_same<stepper_tag, ::pressio::ode::implicitmethods::CrankNicolson>::value and
  ::pressio::containers::predicates::is_vector_wrapper_eigen<residual_type>::value
  >
time_discrete_residual(const fom_states_manager_t & fomStatesMngr,
                       const fom_velocities_manager_t & fomVelocitiesMngr,
		       residual_type & R,
		       const scalar_type & dt,
		       const hyp_ind_t & hypIndices)
{
  using cnst = ::pressio::ode::constants::cranknicolson<scalar_type>;
  constexpr auto cnp1  = cnst::c_np1_;
  constexpr auto cn    = cnst::c_n_;
  constexpr auto cfn   = cnst::c_fn_;
  const auto cfnDt   = cfn*dt;
  const auto cfnp1Dt = cfn*dt;

  const auto & y_np1  = fomStatesMngr(::pressio::ode::nPlusOne());
  const auto & y_n    = fomStatesMngr(::pressio::ode::n());

  const auto & f_np1  = fomVelocitiesMngr(::pressio::ode::nPlusOne());
  const auto & f_n    = fomVelocitiesMngr(::pressio::ode::n());

  assert(R.extent(0) == hypIndices.extent(0));
  for (std::size_t i=0; i<(std::size_t) R.extent(0); ++i)
  {
    const auto yI = hypIndices(i);
    R(i) = cnp1*y_np1(yI) + cn*y_n(yI) + cfnp1Dt*f_np1(i) + cfnDt*f_n(i);
  }
}
#endif

// ----------------------------------------------------------------------
/*
 * for Kokkos and eigen wrappers WITHOUT sample to stencil indices,
 * only works when structures are of same size
 */
// ----------------------------------------------------------------------
template<
  typename stepper_tag,
  typename fom_states_manager_t,
  typename state_type,
  typename scalar_type
  >
void _time_discrete_residual_bdf1_nohr_eig_kokkos(const fom_states_manager_t & fomStatesMngr,
						  state_type & R,
						  scalar_type dt)
{
  const auto & y_np1  = fomStatesMngr.fomStateAt(::pressio::ode::nPlusOne());
  const auto & y_n    = fomStatesMngr.fomStateAt(::pressio::ode::n());

  constexpr auto cnp1 = ::pressio::ode::constants::bdf1<scalar_type>::c_np1_;
  constexpr auto cn   = ::pressio::ode::constants::bdf1<scalar_type>::c_n_;
  const auto cfdt     = ::pressio::ode::constants::bdf1<scalar_type>::c_f_ * dt;
  ::pressio::ops::update(R, cfdt, y_np1, cnp1, y_n, cn);
}

template<
  typename stepper_tag,
  typename fom_states_manager_t,
  typename state_type,
  typename scalar_type
  >
void _time_discrete_residual_bdf2_nohr_eig_kokkos(const fom_states_manager_t & fomStatesMngr,
						  state_type & R,
						  scalar_type dt)
{
  const auto & y_np1 = fomStatesMngr.fomStateAt(::pressio::ode::nPlusOne());
  const auto & y_n   = fomStatesMngr.fomStateAt(::pressio::ode::n());
  const auto & y_nm1 = fomStatesMngr.fomStateAt(::pressio::ode::nMinusOne());

  constexpr auto cnp1 = ::pressio::ode::constants::bdf2<scalar_type>::c_np1_;
  constexpr auto cn   = ::pressio::ode::constants::bdf2<scalar_type>::c_n_;
  constexpr auto cnm1 = ::pressio::ode::constants::bdf2<scalar_type>::c_nm1_;
  const auto cfdt     = ::pressio::ode::constants::bdf2<scalar_type>::c_f_ * dt;

  ::pressio::ops::update(R, cfdt,
			 y_np1, cnp1,
			 y_n, cn,
			 y_nm1, cnm1);
}

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<
  typename stepper_tag,
  typename fom_states_manager_t,
  typename state_type,
  typename scalar_type
  >
::pressio::mpl::enable_if_t<
  std::is_same<stepper_tag, ::pressio::ode::implicitmethods::Euler>::value and
  ::pressio::containers::predicates::is_vector_wrapper_eigen<state_type>::value == true
  >
time_discrete_residual(const fom_states_manager_t & fomStatesMngr,
		       state_type & R,
		       scalar_type dt)
{
  _time_discrete_residual_bdf1_nohr_eig_kokkos<stepper_tag>(fomStatesMngr, R, dt);
}

template<
  typename stepper_tag,
  typename fom_states_manager_t,
  typename state_type,
  typename scalar_type
  >
::pressio::mpl::enable_if_t<
  std::is_same<stepper_tag, ::pressio::ode::implicitmethods::BDF2>::value and
  ::pressio::containers::predicates::is_vector_wrapper_eigen<state_type>::value == true
  >
time_discrete_residual(const fom_states_manager_t & fomStatesMngr,
		       state_type & R,
		       scalar_type dt)
{
  _time_discrete_residual_bdf2_nohr_eig_kokkos<stepper_tag>(fomStatesMngr, R, dt);
}


// Crank-Nicolson
template<
  typename stepper_tag,
  typename fom_states_t,
  typename fom_velocities_t,
  typename state_type,
  typename scalar_type
  >
::pressio::mpl::enable_if_t<
  std::is_same<stepper_tag, ::pressio::ode::implicitmethods::CrankNicolson>::value and
  ::pressio::containers::predicates::is_vector_wrapper_eigen<state_type>::value == true
  >
time_discrete_residual(const fom_states_t & fomStates,
		       const fom_velocities_t & fomVelocities,
		       state_type & R,
		       scalar_type dt)
{
  using cnst = ::pressio::ode::constants::cranknicolson<scalar_type>;
  constexpr auto cnp1  = cnst::c_np1_;
  constexpr auto cn    = cnst::c_n_;
  constexpr auto cfn   = cnst::c_fn_;
  const auto cfnDt   = cfn*dt;
  const auto cfnp1Dt = cfn*dt;

  const auto & y_np1  = fomStates(::pressio::ode::nPlusOne());
  const auto & y_n    = fomStates(::pressio::ode::n());

  ::pressio::ops::update
      (R, y_np1, cnp1, y_n, cn,
       fomVelocities(::pressio::ode::n()), cfnDt,
       fomVelocities(::pressio::ode::nPlusOne()), cfnp1Dt);
}
#endif

#ifdef PRESSIO_ENABLE_TPL_KOKKOS
template<
  typename stepper_tag,
  typename fom_states_manager_t,
  typename state_type,
  typename scalar_type
  >
::pressio::mpl::enable_if_t<
  std::is_same<stepper_tag, ::pressio::ode::implicitmethods::Euler>::value and
  ::pressio::containers::predicates::is_vector_wrapper_kokkos<state_type>::value == true
  >
time_discrete_residual(const fom_states_manager_t & fomStatesMngr,
		       state_type & R,
		       scalar_type dt)
{
  _time_discrete_residual_bdf1_nohr_eig_kokkos<stepper_tag>(fomStatesMngr, R, dt);
}

template<
  typename stepper_tag,
  typename fom_states_manager_t,
  typename state_type,
  typename scalar_type
  >
::pressio::mpl::enable_if_t<
  std::is_same<stepper_tag, ::pressio::ode::implicitmethods::BDF2>::value and
  ::pressio::containers::predicates::is_vector_wrapper_kokkos<state_type>::value == true
  >
time_discrete_residual(const fom_states_manager_t & fomStatesMngr,
		       state_type & R,
		       scalar_type dt)
{
  _time_discrete_residual_bdf2_nohr_eig_kokkos<stepper_tag>(fomStatesMngr, R, dt);
}
#endif


#ifdef PRESSIO_ENABLE_TPL_TRILINOS

/* when we have hyper-reduction and we need to calculate the
 * time-discrete residual, yn (i.e. the current state) and
 * prevStates (i.e. the states at prev steps) might have
 * different entries/maps than R (i.e. the spatial residual)
 * so we need to update R only picking elements corresponding
 * elements in yn and prevStates.
 * In other words, in the functions below we do NOT assume that
 * yn and prevStates have the same distribution/sizes of R but we DO assume
 * that yn and prevStates contain at least all the elements of R (and more).
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
template <typename stepper_tag>
struct time_discrete_single_entry_epetra;

template <>
struct time_discrete_single_entry_epetra<::pressio::ode::implicitmethods::Euler>{
  template <typename T, typename fom_states_manager_t>
    static void evaluate(const T& dt,
			 T & R,
			 int lid,
			 const fom_states_manager_t & fomStatesMngr)
  {
    const auto & y_np1 = fomStatesMngr.fomStateAt(::pressio::ode::nPlusOne());
    const auto & y_n   = fomStatesMngr.fomStateAt(::pressio::ode::n());

    constexpr auto cnp1 = ::pressio::ode::constants::bdf1<T>::c_np1_;
    constexpr auto cn = ::pressio::ode::constants::bdf1<T>::c_n_;
    const auto cfdt	  = ::pressio::ode::constants::bdf1<T>::c_f_ * dt;
    R = cnp1*y_np1(lid) + cn*y_n(lid) + cfdt*R;
  }
};

template <>
struct time_discrete_single_entry_epetra<::pressio::ode::implicitmethods::BDF2>{
  template <typename T, typename fom_states_manager_t>
    static void evaluate(const T& dt,
			 T & R,
			 int lid,
			 const fom_states_manager_t & fomStatesMngr)
  {
    const auto & y_np1 = fomStatesMngr.fomStateAt(::pressio::ode::nPlusOne());
    const auto & y_n   = fomStatesMngr.fomStateAt(::pressio::ode::n());
    const auto & y_nm1 = fomStatesMngr.fomStateAt(::pressio::ode::nMinusOne());

    constexpr auto cnp1   = ::pressio::ode::constants::bdf2<T>::c_np1_;
    constexpr auto cn = ::pressio::ode::constants::bdf2<T>::c_n_;
    constexpr auto cnm1 = ::pressio::ode::constants::bdf2<T>::c_nm1_;
    const auto cfdt	  = ::pressio::ode::constants::bdf2<T>::c_f_ * dt;
    R = cnp1*y_np1(lid) + cn*y_n(lid) + cnm1*y_nm1(lid) + cfdt*R;
  }
};

template<
  typename stepper_tag,
  typename fom_states_manager_t,
  typename state_type,
  typename scalar_type
  >
::pressio::mpl::enable_if_t<
  containers::predicates::is_vector_wrapper_epetra<state_type>::value == true
  >
time_discrete_residual(const fom_states_manager_t & fomStatesMngr,
		       state_type & R,
		       scalar_type dt){
  // On input: R contains the application RHS, i.e. if
  // dudt = f(x,u,...), R contains f(...)

  const auto & y_np1   = fomStatesMngr.fomStateAt(::pressio::ode::nPlusOne());

  // the integral type of the global indices
  using GO_t = typename containers::details::traits<state_type>::global_ordinal_t;

  // map of y_np1 (prevStates has for sure the same map as y_np1)
  const auto & y_map = y_np1.data()->Map();
  // my global elements
  std::vector<GO_t> gIDy( y_np1.extentLocal(0) );
  y_map.MyGlobalElements( gIDy.data() );

  // map of R
  const auto & R_map = R.data()->Map();
  // global elements
  std::vector<GO_t> gIDr( R.extentLocal(0) );
  R_map.MyGlobalElements( gIDr.data() );

  //loop over elements of R
  for (auto i=0; i<R.extentLocal(0); i++){
    // ask the state map what is the local index corresponding
    // to the global index we are handling
    const auto lid = y_map.LID(gIDr[i]);
    // compute the time-discrete entry
    time_discrete_single_entry_epetra<
      stepper_tag
      >::template evaluate<
	scalar_type, fom_states_manager_t>(dt, R(i), lid, fomStatesMngr);
  }
}

/*************************************
	    tpetra
*************************************/
template <typename stepper_tag>
struct time_discrete_single_entry_tpetra;

template <>
struct time_discrete_single_entry_tpetra<::pressio::ode::implicitmethods::Euler>
{
  template <typename T, typename state_type>
    static void evaluate(const T& dt,
			 T & R,
			 int lid,
			 const state_type & ynp1,
			 const state_type & yn)
  {
    constexpr auto cnp1 = ::pressio::ode::constants::bdf1<T>::c_np1_;
    constexpr auto cn   = ::pressio::ode::constants::bdf1<T>::c_n_;
    const auto cfdt	= ::pressio::ode::constants::bdf1<T>::c_f_ * dt;

    R = cnp1*(ynp1.getData())[lid] + cn*(yn.getData())[lid] + cfdt*R;
  }
};

template <>
struct time_discrete_single_entry_tpetra<::pressio::ode::implicitmethods::BDF2>
{
  template <typename T, typename state_type>
    static void evaluate(const T& dt,
			 T & R,
			 int lid,
			 const state_type & ynp1,
			 const state_type & yn,
			 const state_type & ynm1)
  {
    constexpr auto cnp1 = ::pressio::ode::constants::bdf2<T>::c_np1_;
    constexpr auto cn   = ::pressio::ode::constants::bdf2<T>::c_n_;
    constexpr auto cnm1 = ::pressio::ode::constants::bdf2<T>::c_nm1_;
    const auto cfdt	= ::pressio::ode::constants::bdf2<T>::c_f_ * dt;
    R = cnp1*(ynp1.getData())[lid]
      + cn*(yn.getData())[lid]
      + cnm1*(ynm1.getData())[lid]
      + cfdt*R;
  }
};

template<
  typename stepper_tag,
  typename state_type,
  typename scalar_type,
  typename ... Args
  >
::pressio::mpl::enable_if_t<
  containers::predicates::is_vector_tpetra<state_type>::value == true
  >
time_discrete_residual_tpetra_impl(const state_type & currentState,
				   state_type & R,
				   const scalar_type & dt,
				   Args && ... args)
{
  // // On input: R contains the application RHS, i.e. if
  // // dudt = f(x,u,...), R contains f(...)

  // map of currentState (ynm has for sure the same map as currentState)
  const auto y_map = currentState.getMap();
  // my global elements
  const auto gIDy = y_map->getMyGlobalIndices();

  // map of R
  const auto R_map = R.getMap();
  // global elements
  const auto gIDr = R_map->getMyGlobalIndices();

  //loop over elements of R
  for (size_t i=0; i<R.getLocalLength(); i++){
    // ask the state map what is the local index corresponding
    // to the global index we are handling
    const auto lid = y_map->getLocalElement(gIDr[i]);
    // compute the time-discrete entry
    time_discrete_single_entry_tpetra<
      stepper_tag
      >::template evaluate<
	scalar_type, state_type
      >(dt, R.getDataNonConst()[i], lid, currentState, std::forward<Args>(args)...);
  }
}

template<
  typename stepper_tag,
  typename fom_states_manager_t,
  typename state_type,
  typename scalar_type
  >
::pressio::mpl::enable_if_t<
  containers::predicates::is_vector_wrapper_tpetra<state_type>::value == true and
  std::is_same<stepper_tag, ::pressio::ode::implicitmethods::Euler>::value
  >
time_discrete_residual(const fom_states_manager_t & fomStatesMngr,
		       state_type & R,
		       const scalar_type & dt)
{
  const auto & y_np1 = fomStatesMngr.fomStateAt(::pressio::ode::nPlusOne());
  const auto & y_n   = fomStatesMngr.fomStateAt(::pressio::ode::n());
  time_discrete_residual_tpetra_impl<stepper_tag>(*y_np1.data(), *R.data(),
						  dt, *y_n.data());
}

template<
  typename stepper_tag,
  typename fom_states_manager_t,
  typename state_type,
  typename scalar_type
  >
::pressio::mpl::enable_if_t<
  containers::predicates::is_vector_wrapper_tpetra<state_type>::value == true and
  std::is_same<stepper_tag, ::pressio::ode::implicitmethods::BDF2>::value
  >
time_discrete_residual(const fom_states_manager_t & fomStatesMngr,
		       state_type & R,
		       scalar_type dt)
{
  const auto & y_np1 = fomStatesMngr.fomStateAt(::pressio::ode::nPlusOne());
  const auto & y_n   = fomStatesMngr.fomStateAt(::pressio::ode::n());
  const auto & y_nm1 = fomStatesMngr.fomStateAt(::pressio::ode::nMinusOne());

  time_discrete_residual_tpetra_impl<stepper_tag>(*y_np1.data(), *R.data(),
						  dt, *y_n.data(), *y_nm1.data());
}

/*************************************
            tpetra block
*************************************/
template<
  typename stepper_tag,
  typename fom_states_manager_t,
  typename state_type,
  typename scalar_type
  >
::pressio::mpl::enable_if_t<
  containers::predicates::is_vector_wrapper_tpetra_block<state_type>::value and
  std::is_same<stepper_tag, ::pressio::ode::implicitmethods::Euler>::value
  >
time_discrete_residual(const fom_states_manager_t & fomStatesMngr,
		       state_type & R,
		       const scalar_type & dt)
{
  const auto & y_np1 = fomStatesMngr.fomStateAt(::pressio::ode::nPlusOne());
  const auto & y_n   = fomStatesMngr.fomStateAt(::pressio::ode::n());

  auto ynp1_vv = const_cast<state_type &>(y_np1).data()->getVectorView();
  auto yn_vv   = const_cast<state_type &>(y_n).data()->getVectorView();
  auto R_vv    = R.data()->getVectorView();

  time_discrete_residual_tpetra_impl<stepper_tag>(ynp1_vv, R_vv, dt, yn_vv);
}

template<
  typename stepper_tag,
  typename fom_states_manager_t,
  typename state_type,
  typename scalar_type
  >
::pressio::mpl::enable_if_t<
  containers::predicates::is_vector_wrapper_tpetra_block<state_type>::value and
  std::is_same<stepper_tag, ::pressio::ode::implicitmethods::BDF2>::value
  >
time_discrete_residual(const fom_states_manager_t & fomStatesMngr,
		       state_type & R,
		       scalar_type dt)
{
  const auto & y_np1 = fomStatesMngr.fomStateAt(::pressio::ode::nPlusOne());
  const auto & y_n   = fomStatesMngr.fomStateAt(::pressio::ode::n());
  const auto & y_nm1 = fomStatesMngr.fomStateAt(::pressio::ode::nMinusOne());

  auto ynp1_vv = const_cast<state_type &>(y_np1).data()->getVectorView();
  auto yn_vv   = const_cast<state_type &>(y_n).data()->getVectorView();
  auto ynm1_vv = const_cast<state_type &>(y_nm1).data()->getVectorView();
  auto R_vv = R.data()->getVectorView();
  time_discrete_residual_tpetra_impl<stepper_tag>(ynp1_vv, R_vv, dt, yn_vv, ynm1_vv);
}

#endif

}}}}}
#endif  // ROM_LSPG_IMPL_UNSTEADY_CONTINUOUS_TIME_API_DISCRETE_TIME_FUNCTIONS_ROM_LSPG_TIME_DISCRETE_RESIDUAL_HPP_
