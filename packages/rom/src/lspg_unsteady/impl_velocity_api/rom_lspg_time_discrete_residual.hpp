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

#ifndef ROM_LSPG_TIME_DISCRETE_RESIDUAL_HPP_
#define ROM_LSPG_TIME_DISCRETE_RESIDUAL_HPP_

namespace pressio{ namespace rom{ namespace lspg{ namespace unsteady{ namespace impl{

/* enable when we have:
 * regular c++, BDF1 and user-defined ops
*/
template<
  typename stepper_tag,
  typename fom_states_cont_t,
  typename state_type,
  typename scalar_type,
  typename ud_ops,
  ::pressio::mpl::enable_if_t<
    std::is_same<stepper_tag, ::pressio::ode::implicitmethods::Euler>::value
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
    and !::pressio::containers::meta::is_vector_wrapper_pybind<state_type>::value
    and mpl::not_same< ud_ops, pybind11::object>::value
#endif
   > * = nullptr
  >
void time_discrete_residual(const fom_states_cont_t & fomStates,
			    state_type & R,
			    const scalar_type & dt,
			    const ud_ops * udOps){

  const auto & fomStateAt_n   = fomStates.getCRefToCurrentFomState();
  const auto & fomStateAt_nm1 = fomStates.getCRefToFomStatePrevStep();
  udOps->time_discrete_euler(*R.data(), *fomStateAt_n.data(), *fomStateAt_nm1.data(), dt);
}


#ifdef PRESSIO_ENABLE_TPL_PYBIND11
/*
 * for python binddings, enable when we have BDF1 and we do the computation
*/
template<
  typename stepper_tag,
  typename fom_states_cont_t,
  typename residual_type,
  typename scalar_type,
  ::pressio::mpl::enable_if_t<
    std::is_same<stepper_tag, ::pressio::ode::implicitmethods::Euler>::value and
    ::pressio::containers::meta::is_vector_wrapper_pybind<residual_type>::value
    > * = nullptr
  >
void time_discrete_residual(const fom_states_cont_t & fomStates,
			    residual_type & R,
			    const scalar_type & dt)
{
  auto & fomStateAt_n   = fomStates.getCRefToCurrentFomState();
  auto & fomStateAt_nm1 = fomStates.getCRefToFomStatePrevStep();

  constexpr auto cn   = ::pressio::ode::constants::bdf1<scalar_type>::c_n_;
  constexpr auto cnm1 = ::pressio::ode::constants::bdf1<scalar_type>::c_nm1_;
  const auto cf	      = ::pressio::ode::constants::bdf1<scalar_type>::c_f_ * dt;

  // //R = y_n - y_nm1 - dt * R;
  ::pressio::ops::do_update(R, cf, fomStateAt_n, cn, fomStateAt_nm1, cnm1);
}
#endif



// ----------------------------------------------------------------------
/*
 * for Kokkos and eigen wrappers (this is not suitable for sample mesh)
 * only works when structures are of same size
*/
// ----------------------------------------------------------------------
template<
  typename stepper_tag,
  typename fom_states_cont_t,
  typename state_type,
  typename scalar_type,
  ::pressio::mpl::enable_if_t<
    std::is_same<stepper_tag, ::pressio::ode::implicitmethods::Euler>::value and
    (containers::meta::is_vector_wrapper_eigen<state_type>::value == true
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
     or containers::meta::is_vector_wrapper_kokkos<state_type>::value == true
#endif
     )
    > * = nullptr
  >
void time_discrete_residual(const fom_states_cont_t & fomStates,
			    state_type & R,
			    scalar_type dt){

  const auto & fomStateAt_n   = fomStates.getCRefToCurrentFomState();
  const auto & fomStateAt_nm1 = fomStates.getCRefToFomStatePrevStep();

  assert( R.extent(0) == fomStateAt_n.extent(0) );
  assert( fomStateAt_n.extent(0) == fomStateAt_nm1.extent(0) );

  constexpr auto cn   = ::pressio::ode::constants::bdf1<scalar_type>::c_n_;
  constexpr auto cnm1 = ::pressio::ode::constants::bdf1<scalar_type>::c_nm1_;
  const auto cf	  = ::pressio::ode::constants::bdf1<scalar_type>::c_f_ * dt;

  //R = y_n - y_nm1 - dt * R;
  ::pressio::ops::do_update(R, cf, fomStateAt_n, cn, fomStateAt_nm1, cnm1);
}


template<
  typename stepper_tag,
  typename fom_states_cont_t,
  typename state_type,
  typename scalar_type,
  ::pressio::mpl::enable_if_t<
    std::is_same<stepper_tag, ::pressio::ode::implicitmethods::BDF2>::value and
    (containers::meta::is_vector_wrapper_eigen<state_type>::value == true
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
     or containers::meta::is_vector_wrapper_kokkos<state_type>::value == true
#endif
     )
    > * = nullptr
  >
void time_discrete_residual(const fom_states_cont_t & fomStates,
			    state_type		& R,
			    const scalar_type	& dt){

  const auto & fomStateAt_n   = fomStates.getCRefToCurrentFomState();
  const auto & fomStateAt_nm1   = fomStates.getCRefToFomStatePrevStep();
  const auto & fomStateAt_nm2   = fomStates.getCRefToFomStatePrevPrevStep();

  assert( R.extent(0) == fomStateAt_n.extent(0) );
  assert( fomStateAt_n.extent(0) == fomStateAt_nm1.extent(0) );
  assert( fomStateAt_nm1.extent(0) == fomStateAt_nm2.extent(0));

  constexpr auto cn   = ::pressio::ode::constants::bdf2<scalar_type>::c_n_;
  constexpr auto cnm1 = ::pressio::ode::constants::bdf2<scalar_type>::c_nm1_;
  constexpr auto cnm2 = ::pressio::ode::constants::bdf2<scalar_type>::c_nm2_;
  const auto cf	  = ::pressio::ode::constants::bdf2<scalar_type>::c_f_ * dt;

  auto & y_nm1 = fomStates.getCRefToFomStatePrevStep();
  auto & y_nm2 = fomStates.getCRefToFomStatePrevPrevStep();
  // compute: R = y_n - 4/3 * y_n-1 + 1/3 * y_n-2 - 2/3 * dt * f(y_n, t_n)
  ::pressio::ops::do_update(R, cf, fomStateAt_n, cn, y_nm1, cnm1, y_nm2, cnm2);
}





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
  template <typename T, typename fom_states_cont_t>
    static void evaluate(const T& dt,
			 T & R,
			 int lid,
			 const fom_states_cont_t & odeStates)
  {
    const auto & y_n   = odeStates.getCRefToCurrentFomState();
    const auto & y_nm1 = odeStates.getCRefToFomStatePrevStep();

    constexpr auto cn   = ::pressio::ode::constants::bdf1<T>::c_n_;
    constexpr auto cnm1 = ::pressio::ode::constants::bdf1<T>::c_nm1_;
    const auto cf	  = ::pressio::ode::constants::bdf1<T>::c_f_ * dt;
    R = cn*y_n[lid] + cnm1*y_nm1[lid] + cf*R;
  }
};

template <>
struct time_discrete_single_entry_epetra<::pressio::ode::implicitmethods::BDF2>{
  template <typename T, typename fom_states_cont_t>
    static void evaluate(const T& dt,
			 T & R,
			 int lid,
			 const fom_states_cont_t & odeStates)
  {
    const auto & y_n   = odeStates.getCRefToCurrentFomState();
    const auto & y_nm1 = odeStates.getCRefToFomStatePrevStep();
    const auto & y_nm2 = odeStates.getCRefToFomStatePrevPrevStep();

    constexpr auto cn   = ::pressio::ode::constants::bdf2<T>::c_n_;
    constexpr auto cnm1 = ::pressio::ode::constants::bdf2<T>::c_nm1_;
    constexpr auto cnm2 = ::pressio::ode::constants::bdf2<T>::c_nm2_;
    const auto cf	  = ::pressio::ode::constants::bdf2<T>::c_f_ * dt;
    R = cn*y_n[lid] + cnm1*y_nm1[lid] + cnm2*y_nm2[lid] + cf*R;
  }
};


template<
  typename stepper_tag,
  typename fom_states_cont_t,
  typename state_type,
  typename scalar_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_vector_wrapper_epetra<state_type>::value == true
    > * = nullptr
  >
void time_discrete_residual(const fom_states_cont_t & fomStates,
			    state_type & R,
			    scalar_type dt){
  // On input: R contains the application RHS, i.e. if
  // dudt = f(x,u,...), R contains f(...)

  const auto & y_n   = fomStates.getCRefToCurrentFomState();

  // the integral type of the global indices
  using GO_t = typename containers::details::traits<state_type>::global_ordinal_t;

  // get map of y_n (prevStates has for sure the same map as y_n)
  const auto & y_map = y_n.data()->Map();
  // get my global elements
  std::vector<GO_t> gIDy( y_n.extentLocal(0) );
  y_map.MyGlobalElements( gIDy.data() );

  // get map of R
  const auto & R_map = R.data()->Map();
  // get global elements
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
      scalar_type, fom_states_cont_t>(dt, R[i], lid, fomStates);
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
		       const state_type & yn,
		       const state_type & ynm1){

    constexpr auto cn   = ::pressio::ode::constants::bdf1<T>::c_n_;
    constexpr auto cnm1 = ::pressio::ode::constants::bdf1<T>::c_nm1_;
    const auto cf	= ::pressio::ode::constants::bdf1<T>::c_f_ * dt;

    R = cn*(yn.getData())[lid] + cnm1*(ynm1.getData())[lid] + cf*R;
  }
};

template <>
struct time_discrete_single_entry_tpetra<::pressio::ode::implicitmethods::BDF2>
{
  template <typename T, typename state_type>
    static void evaluate(const T& dt,
			 T & R,
			 int lid,
			 const state_type& yn,
			 const state_type & ynm1,
			 const state_type & ynm2)
  {
    constexpr auto cn   = ::pressio::ode::constants::bdf2<T>::c_n_;
    constexpr auto cnm1 = ::pressio::ode::constants::bdf2<T>::c_nm1_;
    constexpr auto cnm2 = ::pressio::ode::constants::bdf2<T>::c_nm2_;
    const auto cf	  = ::pressio::ode::constants::bdf2<T>::c_f_ * dt;
    R = cn*(yn.getData())[lid]
      + cnm1*(ynm1.getData())[lid]
      + cnm2*(ynm2.getData())[lid]
      + cf*R;
  }
};



template<
  typename stepper_tag,
  typename state_type,
  typename scalar_type,
  typename ... Args,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_vector_tpetra<state_type>::value == true
    > * = nullptr
  >
void time_discrete_residual_tpetra_impl(const state_type & currentState,
					state_type & R,
					const scalar_type & dt,
					Args && ... args)
{
  // // On input: R contains the application RHS, i.e. if
  // // dudt = f(x,u,...), R contains f(...)

  // get map of currentState (ynm has for sure the same map as currentState)
  const auto y_map = currentState.getMap();
  // get my global elements
  const auto gIDy = y_map->getMyGlobalIndices();

  // get map of R
  const auto R_map = R.getMap();
  // get global elements
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
  typename fom_states_cont_t,
  typename state_type,
  typename scalar_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_vector_wrapper_tpetra<state_type>::value == true and
    std::is_same<stepper_tag, ::pressio::ode::implicitmethods::Euler>::value
    > * = nullptr
  >
void time_discrete_residual(const fom_states_cont_t & odeStates,
			    state_type & R,
			    const scalar_type & dt){

  const auto & y_n   = odeStates.getCRefToCurrentFomState();
  const auto & y_nm1 = odeStates.getCRefToFomStatePrevStep();
  time_discrete_residual_tpetra_impl<stepper_tag>(*y_n.data(), *R.data(), dt, *y_nm1.data());

}

template<
  typename stepper_tag,
  typename fom_states_cont_t,
  typename state_type,
  typename scalar_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_vector_wrapper_tpetra<state_type>::value == true and
    std::is_same<stepper_tag, ::pressio::ode::implicitmethods::BDF2>::value
    > * = nullptr
  >
void time_discrete_residual(const fom_states_cont_t & fomStates,
			    state_type & R,
			    scalar_type dt){

  const auto & y_n   = fomStates.getCRefToCurrentFomState();
  const auto & y_nm1 = fomStates.getCRefToFomStatePrevStep();
  const auto & y_nm2 = fomStates.getCRefToFomStatePrevPrevStep();

  time_discrete_residual_tpetra_impl<stepper_tag>(*y_n.data(), *R.data(), dt,
						     *y_nm1.data(),
						     *y_nm2.data());
}


/*************************************
            tpetra block
*************************************/
template<
  typename stepper_tag,
  typename fom_states_cont_t,
  typename state_type,
  typename scalar_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_vector_wrapper_tpetra_block<state_type>::value and
    std::is_same<stepper_tag, ::pressio::ode::implicitmethods::Euler>::value
    > * = nullptr
  >
void time_discrete_residual(const fom_states_cont_t & fomStates,
			    state_type & R,
			    const scalar_type & dt)
{
  const auto & y_n   = fomStates.getCRefToCurrentFomState();
  const auto & y_nm1 = fomStates.getCRefToFomStatePrevStep();

  auto yn_vv   = const_cast<state_type &>(y_n).data()->getVectorView();
  auto ynm0_vv = const_cast<state_type &>(y_nm1).data()->getVectorView();
  auto R_vv    = R.data()->getVectorView();

  time_discrete_residual_tpetra_impl<stepper_tag>(yn_vv, R_vv, dt, ynm0_vv);
}

template<
  typename stepper_tag,
  typename fom_states_cont_t,
  typename state_type,
  typename scalar_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_vector_wrapper_tpetra_block<state_type>::value and
    std::is_same<stepper_tag, ::pressio::ode::implicitmethods::BDF2>::value
    > * = nullptr
  >
void time_discrete_residual(const fom_states_cont_t & fomStates,
			    state_type & R,
			    scalar_type dt){

  const auto & y_n   = fomStates.getCRefToCurrentFomState();
  const auto & y_nm1 = fomStates.getCRefToFomStatePrevStep();
  const auto & y_nm2 = fomStates.getCRefToFomStatePrevPrevStep();

  auto yn_vv   = const_cast<state_type &>(y_n).data()->getVectorView();
  auto ynm1_vv = const_cast<state_type &>(y_nm1).data()->getVectorView();
  auto ynm2_vv = const_cast<state_type &>(y_nm2).data()->getVectorView();
  auto R_vv = R.data()->getVectorView();
  time_discrete_residual_tpetra_impl<stepper_tag>(yn_vv, R_vv, dt, ynm1_vv, ynm2_vv);
}

#endif

}}}}}//end namespace pressio::rom::lspg::unstedy::impl
#endif
