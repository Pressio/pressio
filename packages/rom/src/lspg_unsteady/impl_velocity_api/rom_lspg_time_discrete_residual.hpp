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

#include "../../../../ode/src/ode_ConfigDefs.hpp"
#include "../../../../ode/src/implicit/ode_implicit_constants.hpp"

namespace pressio{ namespace rom{ namespace impl{

template<
  ::pressio::ode::ImplicitEnum method,
  int n,
  typename state_type,
  typename scalar_type,
  typename ud_ops,
  ::pressio::mpl::enable_if_t<
    method == ::pressio::ode::ImplicitEnum::Euler
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
    and mpl::not_same< ud_ops, pybind11::object>::value
#endif
   > * = nullptr
  >
void time_discrete_residual(const state_type & yn,
			    const std::array<state_type,n> & ynm,
			    state_type & R,
			    scalar_type dt,
			    const ud_ops * udOps){

  udOps->time_discrete_euler(*R.data(), *yn.data(), *ynm[0].data(), dt);
}

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
template<
  ::pressio::ode::ImplicitEnum method,
  int n,
  typename state_type,
  typename scalar_type,
  typename ud_ops,
  ::pressio::mpl::enable_if_t<
    method == ::pressio::ode::ImplicitEnum::Euler and
    mpl::is_same< ud_ops, pybind11::object>::value
    > * = nullptr
  >
void time_discrete_residual(const state_type & yn,
			    const std::array<state_type,n> & ynm,
			    state_type & R,
			    scalar_type dt,
			    const ud_ops & udOps){

  //TODO: this function name will need to be changed to something better
  udOps.attr("time_discrete_euler")(R, yn, ynm[0], dt);
}
#endif



/*
 * for Kokkos and eigen wrappers (this is not suitable for sample mesh)
 * only works when structures are of same size
*/
template<
  ::pressio::ode::ImplicitEnum method,
  int n,
  typename state_type,
  typename scalar_type,
  ::pressio::mpl::enable_if_t<
    method == ::pressio::ode::ImplicitEnum::Euler and
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
    (containers::meta::is_vector_wrapper_kokkos<state_type>::value == true
     or
#endif
     containers::meta::is_vector_wrapper_eigen<state_type>::value == true
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
     )
#endif
    > * = nullptr
  >
void time_discrete_residual(const state_type & yn,
			    const std::array<state_type,n> & ynm,
			    state_type & R,
			    scalar_type dt){

  assert( R.size() == yn.size() );
  assert( yn.size() == ynm[0].size() );
  constexpr auto one = ::pressio::utils::constants::one<scalar_type>();
  constexpr auto negOne = ::pressio::utils::constants::negOne<scalar_type>();
  const auto negDt = negOne*dt;
  //R = yn - ynm[0] - dt * R;
  ::pressio::containers::ops::do_update(R, negDt, yn, one, ynm[0], negOne);
}


template<
  ::pressio::ode::ImplicitEnum method,
  int n,
  typename state_type,
  typename scalar_type,
  ::pressio::mpl::enable_if_t<
    method == ::pressio::ode::ImplicitEnum::BDF2 and
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
    (containers::meta::is_vector_wrapper_kokkos<state_type>::value == true
     or
#endif
     containers::meta::is_vector_wrapper_eigen<state_type>::value == true
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
     )
#endif
    > * = nullptr
  >
void time_discrete_residual(const state_type & yn,
			    const std::array<state_type,n> & ynm,
			    state_type & R,
			    scalar_type dt){
  assert( R.size() == yn.size() );
  assert( yn.size() == ynm[0].size() );
  assert( ynm[0].size() == ynm[1].size());

  using namespace ::pressio::ode::coeffs;

  constexpr auto one = ::pressio::utils::constants::one<scalar_type>();
  constexpr auto negOne = ::pressio::utils::constants::negOne<scalar_type>();
  const auto a = negOne * bdf2<scalar_type>::c3_*dt;
  const auto b = bdf2<scalar_type>::c2_;
  const auto c = negOne * bdf2<scalar_type>::c1_;

  ::pressio::containers::ops::do_update(R, a, yn, one, ynm[0], b, ynm[1], c);
  //   R = yn
  //     - bdf2<scalar_type>::c1_*ynm[1]
  //     + bdf2<scalar_type>::c2_*ynm[0]
  //     - bdf2<scalar_type>::c3_*dt*R;
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

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
template <::pressio::ode::ImplicitEnum odeMethod>
struct time_discrete_single_entry_epetra;

template <>
struct time_discrete_single_entry_epetra<::pressio::ode::ImplicitEnum::Euler>{
  template <typename T, typename state_type, int n>
    static void evaluate(const T& dt,
			 T & R,
			 int lid,
			 const state_type& yn,
			 const std::array<state_type,n> & ynm){
    R = yn[lid] - ynm[0][lid] - dt * R;
  }
};

template <>
struct time_discrete_single_entry_epetra<::pressio::ode::ImplicitEnum::BDF2>{
  template <typename T, typename state_type, int n>
    static void evaluate(const T& dt, T & R,
			 int lid,
			 const state_type& yn,
			 const std::array<state_type,n> & ynm){
    using namespace ::pressio::ode::coeffs;

    R = yn[lid]
      - bdf2<T>::c1_*ynm[1][lid]
      + bdf2<T>::c2_*ynm[0][lid]
      - bdf2<T>::c3_*dt*R;
  }
};


template<
  ::pressio::ode::ImplicitEnum odeMethod,
  int numStates,
  typename state_type,
  typename scalar_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_vector_wrapper_epetra<state_type>::value == true
    > * = nullptr
  >
void time_discrete_residual(const state_type & yn,
			    const std::array<state_type,numStates> & ynm,
			    state_type & R,
			    scalar_type dt){
  // On input: R contains the application RHS, i.e. if
  // dudt = f(x,u,...), R contains f(...)

  // the integral type of the global indices
  using GO_t = typename containers::details::traits<state_type>::global_ordinal_t;

  // get map of yn (ynm has for sure the same map as yn)
  const auto & y_map = yn.getDataMap();
  // get my global elements
  std::vector<GO_t> gIDy( yn.localSize() );
  y_map.MyGlobalElements( gIDy.data() );

  // get map of R
  const auto & R_map = R.getDataMap();
  // get global elements
  std::vector<GO_t> gIDr( R.localSize() );
  R_map.MyGlobalElements( gIDr.data() );

  //loop over elements of R
  for (auto i=0; i<R.localSize(); i++){
    // ask the state map what is the local index corresponding
    // to the global index we are handling
    auto lid = y_map.LID(gIDr[i]);
    // compute the time-discrete entry
    time_discrete_single_entry_epetra<
      odeMethod
      >::template evaluate<
      scalar_type, state_type, numStates>(dt, R[i], lid, yn, ynm);
  }
}




/*************************************
	             tpetra
*************************************/
template <::pressio::ode::ImplicitEnum odeMethod>
struct time_discrete_single_entry_tpetra;

template <>
struct time_discrete_single_entry_tpetra<::pressio::ode::ImplicitEnum::Euler>{
  template <typename T, typename state_type, int n>
    static void evaluate(const T& dt, T & R,
			 int lid,
			 const state_type& yn,
			 const std::array<state_type,n> & ynm){
    //R[i] = yn[lid]- ynm[0][lid] - dt*R[i];
    R = (yn.getData())[lid]
      - (ynm[0].getData())[lid] - dt * R;
  }

  template <typename T, typename state_type>
    static void evaluate(const T& dt, T & R,
			 int lid,
			 const state_type& yn,
			 const state_type & ynm0){
    R = (yn.getData())[lid]
      - (ynm0.getData())[lid] - dt * R;
  }
};

template <>
struct time_discrete_single_entry_tpetra<::pressio::ode::ImplicitEnum::BDF2>{
  template <typename T, typename state_type, int n>
    static void evaluate(const T& dt,
			 T & R,
			 int lid,
			 const state_type & yn,
			 const std::array<state_type,n> & ynm){
    using namespace ::pressio::ode::coeffs;

    R = (yn.getData())[lid]
      - bdf2<T>::c1_*(ynm[1].getData())[lid]
      + bdf2<T>::c2_*(ynm[0].getData())[lid]
      - bdf2<T>::c3_*dt*R;
  }

  template <typename T, typename state_type>
    static void evaluate(const T& dt,
			 T & R,
			 int lid,
			 const state_type& yn,
			 const state_type & ynm0,
			 const state_type & ynm1){
    using namespace ::pressio::ode::coeffs;
    R = (yn.getData())[lid]
      - bdf2<T>::c1_*(ynm1.getData())[lid]
      + bdf2<T>::c2_*(ynm0.getData())[lid]
      - bdf2<T>::c3_*dt*R;
  }
};


template<
  ::pressio::ode::ImplicitEnum odeMethod,
  typename state_type,
  typename scalar_type,
  typename ... Args,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_vector_tpetra<state_type>::value == true
    > * = nullptr
  >
void time_discrete_residual_tpetra_impl(const state_type & yn,
			    state_type & R,
			    scalar_type dt,
			    Args && ... args)/*
			    const state_type & ynm0,
			    const state_type & ynm1,)*/
{
  // // On input: R contains the application RHS, i.e. if
  // // dudt = f(x,u,...), R contains f(...)

  // get map of yn (ynm has for sure the same map as yn)
  const auto y_map = yn.getMap();
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
      odeMethod
      >::template evaluate<
      scalar_type, state_type
      >(dt, R.getDataNonConst()[i], lid, yn,
	std::forward<Args>(args)...);//ynm0, ynm1);
  }
}


template<
  ::pressio::ode::ImplicitEnum odeMethod,
  int n,
  typename state_type,
  typename scalar_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_vector_wrapper_tpetra<state_type>::value == true and
    odeMethod == ::pressio::ode::ImplicitEnum::Euler
    > * = nullptr
  >
void time_discrete_residual(const state_type & yn,
			    const std::array<state_type,n> & ynm,
			    state_type & R,
			    scalar_type dt){

  time_discrete_residual_tpetra_impl<odeMethod>
    (*yn.data(),
     *R.data(), dt,
     *ynm[0].data());

}

template<
  ::pressio::ode::ImplicitEnum odeMethod,
  int n,
  typename state_type,
  typename scalar_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_vector_wrapper_tpetra<state_type>::value == true and
    odeMethod == ::pressio::ode::ImplicitEnum::BDF2
    > * = nullptr
  >
void time_discrete_residual(const state_type & yn,
			    const std::array<state_type,n> & ynm,
			    state_type & R,
			    scalar_type dt){

  time_discrete_residual_tpetra_impl<odeMethod>
    (*yn.data(),
     *R.data(), dt,
     *ynm[0].data(),
     *ynm[1].data());
}


/*************************************
            tpetra block
*************************************/
template<
  ::pressio::ode::ImplicitEnum odeMethod,
  int n,
  typename state_type,
  typename scalar_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_vector_wrapper_tpetra_block<state_type>::value and
    odeMethod == ::pressio::ode::ImplicitEnum::Euler
    > * = nullptr
  >
void time_discrete_residual(const state_type & yn,
			    const std::array<state_type,n> & ynm,
			    state_type & R,
			    scalar_type dt){
  auto yn_vv = const_cast<state_type &>(yn).data()->getVectorView();
  auto ynm0_vv = const_cast<state_type &>(ynm[0]).data()->getVectorView();
  auto R_vv = R.data()->getVectorView();
  time_discrete_residual_tpetra_impl<odeMethod>(
      yn_vv, R_vv, dt, ynm0_vv);
}

template<
  ::pressio::ode::ImplicitEnum odeMethod,
  int n,
  typename state_type,
  typename scalar_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_vector_wrapper_tpetra_block<state_type>::value and
    odeMethod == ::pressio::ode::ImplicitEnum::BDF2
    > * = nullptr
  >
void time_discrete_residual(const state_type & yn,
			    const std::array<state_type,n> & ynm,
			    state_type & R,
			    scalar_type dt){
  auto yn_vv = const_cast<state_type &>(yn).data()->getVectorView();
  auto ynm0_vv = const_cast<state_type &>(ynm[0]).data()->getVectorView();
  auto ynm1_vv = const_cast<state_type &>(ynm[1]).data()->getVectorView();
  auto R_vv = R.data()->getVectorView();
  time_discrete_residual_tpetra_impl<odeMethod>(
    yn_vv, R_vv, dt, ynm0_vv, ynm1_vv);
}

#endif

}}}//end namespace pressio::rom::impl
#endif
