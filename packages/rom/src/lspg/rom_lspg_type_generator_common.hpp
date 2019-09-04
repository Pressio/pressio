/*
//@HEADER
// ************************************************************************
//
// rom_lspg_type_generator_common.hpp
//                     		      Pressio 
// Copyright 2019 National Technology & Engineering Solutions of Sandia,LLC 
//							      (NTESS)
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

#ifndef ROM_LSPG_TYPE_GENERATOR_COMMON_HPP_
#define ROM_LSPG_TYPE_GENERATOR_COMMON_HPP_

#include "../rom_ConfigDefs.hpp"
#include "../rom_fwd.hpp"
#include "../rom_data_fom_rhs.hpp"
#include "../rom_data_fom_states.hpp"
#include "../rom_reconstructor_fom_state.hpp"
#include "../policies/rom_evaluate_fom_velocity_steady_policy.hpp"
#include "../policies/rom_evaluate_fom_velocity_unsteady_policy.hpp"
#include "../policies/rom_apply_fom_jacobian_steady_policy.hpp"
#include "../policies/rom_apply_fom_jacobian_unsteady_policy.hpp"
#include "../../../ode/src/ode_fwd.hpp"

namespace pressio{ namespace rom{

template <ode::ImplicitEnum odeName>
struct statesStorageCapacityHelper{
  static constexpr int maxAuxStates_ = 1;
};

template <>
struct statesStorageCapacityHelper<ode::ImplicitEnum::Euler>{
  static constexpr int maxAuxStates_ = 1;
};

template <>
struct statesStorageCapacityHelper<ode::ImplicitEnum::BDF2>{
  static constexpr int maxAuxStates_ = 2;
};
//-------------------------------------------------------


template <ode::ImplicitEnum odeName, typename ... Rest>
struct auxStepperHelper{
  using type = void;
};

template <
typename ... Rest
>
struct auxStepperHelper<
  ode::ImplicitEnum::BDF2, Rest...
  >{
  using type = ode::ImplicitStepper<
    ode::ImplicitEnum::Euler, Rest...>;
};


//-------------------------------------------------------

template <
  typename fom_type,
  typename decoder_type,
  typename lspg_state_type,
  ode::ImplicitEnum odeName = ode::ImplicitEnum::Undefined,
  typename ud_ops = void,
  typename enable = void
  >
struct LSPGCommonTypes;


//-------------------------------------------------------
// partial specialize for when we are native C++
//-------------------------------------------------------
template <
  typename fom_type,
  typename decoder_type,
  typename lspg_state_type,
  ode::ImplicitEnum odeName,
  typename ud_ops
  >
struct LSPGCommonTypes<
  fom_type, decoder_type, lspg_state_type,
  odeName, ud_ops,
  mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper<lspg_state_type>::value
#ifdef HAVE_PYBIND11
    and mpl::not_same<fom_type, pybind11::object>::value
#endif
    >
  >
{
  LSPGCommonTypes() = default;
  ~LSPGCommonTypes() = default;

  // these are native types of the full-order model (fom)
  using fom_t			= fom_type;
  using scalar_t		= typename fom_t::scalar_type;
  using fom_native_state_t	= typename fom_t::state_type;
  using fom_native_velocity_t	= typename fom_t::velocity_type;

  // fom wrapper types
  using fom_state_t	= ::pressio::containers::Vector<fom_native_state_t>;
  using fom_velocity_t	= ::pressio::containers::Vector<fom_native_velocity_t>;

  // rom state type (passed in)
  using lspg_state_t		= lspg_state_type;

  // for LSPG, the rom residual type = containers::wrapper of application rhs
  // i.e. the wrapped fom rhs type
  using lspg_residual_t		= fom_velocity_t;

  // decoder types (passed in)
  using decoder_t		= decoder_type;
  using decoder_jac_t		= typename decoder_t::jacobian_t;

  // fom state reconstructor type
  using fom_state_reconstr_t	= FomStateReconstructor<
    fom_state_t, decoder_t>;

  // max num of states needed for time integration.
  // this is deduced based on the integrator, i.e. odeName
  static constexpr auto maxAuxStates =
    statesStorageCapacityHelper<odeName>::maxAuxStates_;

  // class type holding fom states data
  using fom_states_data = ::pressio::rom::FomStatesData<
	fom_state_t, maxAuxStates, fom_state_reconstr_t>;

  // class type holding fom rhs data
  using fom_velocity_data = ::pressio::rom::FomRhsData<fom_velocity_t>;

  // if we have a non-trivial user-defined ops
  using ud_ops_t = ud_ops;
};



#ifdef HAVE_PYBIND11
//-------------------------------------------------------
// partial specialize for pybind11
//-------------------------------------------------------
template <
  typename fom_type,
  typename decoder_type,
  typename lspg_state_type,
  ode::ImplicitEnum odeName,
  typename ud_ops
  >
struct LSPGCommonTypes<
  fom_type, decoder_type, lspg_state_type,
  odeName, ud_ops,
  mpl::enable_if_t<
    ::pressio::containers::meta::is_cstyle_array_pybind11<lspg_state_type>::value and
    mpl::is_same<fom_type, pybind11::object>::value
    >
  >
{
  // in this case there is no difference between types because
  // they all are pybind11::array_t so basically wrappers of numpy arrays
  // Since this is used to interface to python, EVERYTHING is done using numpy arrays

  LSPGCommonTypes() = default;
  ~LSPGCommonTypes() = default;

  // these are native types of the full-order model (fom)
  using fom_t			= fom_type;
  using scalar_t		= typename decoder_type::scalar_t;
  using fom_native_state_t	= lspg_state_type;
  using fom_native_velocity_t	= lspg_state_type;

  // fom types
  using fom_state_t	= lspg_state_type;
  using fom_velocity_t	= lspg_state_type;

  // rom state type (passed in)
  using lspg_state_t		= lspg_state_type;

  // for LSPG, the rom residual
  using lspg_residual_t		= fom_velocity_t;

  // decoder types (passed in)
  using decoder_t		= decoder_type;
  using decoder_jac_t		= typename decoder_t::jacobian_t;

  // fom state reconstructor type
  using fom_state_reconstr_t	= FomStateReconstructor<
    fom_state_t, decoder_t>;

  // max num of states needed for time integration.
  // this is deduced based on the integrator, i.e. odeName
  static constexpr auto maxAuxStates =
    statesStorageCapacityHelper<odeName>::maxAuxStates_;

  // class type holding fom states data
  using fom_states_data = ::pressio::rom::FomStatesData<
	fom_state_t, maxAuxStates, fom_state_reconstr_t>;

  // class type holding fom rhs data
  using fom_velocity_data = ::pressio::rom::FomRhsData<fom_velocity_t>;

  // if we have a non-trivial user-defined ops
  using ud_ops_t = ud_ops;
};
#endif


}}//end  namespace pressio::rom
#endif
