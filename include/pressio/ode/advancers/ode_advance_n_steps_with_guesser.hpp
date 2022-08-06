/*
//@HEADER
// ************************************************************************
//
// ode_advance_n_steps.hpp
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

#ifndef ODE_ADVANCERS_ODE_ADVANCE_N_STEPS_WITH_GUESSER_HPP_
#define ODE_ADVANCERS_ODE_ADVANCE_N_STEPS_WITH_GUESSER_HPP_

#include "./impl/ode_advance_noop_observer.hpp"
#include "./impl/ode_advance_n_steps.hpp"

namespace pressio{ namespace ode{

// ---------------------------------
// overload set for stepper
// ---------------------------------

template<
  class StepperType, class StateType, class GuesserType, class IndVarType
  >
mpl::enable_if_t<
  Steppable<StepperType>::value
  && StateGuesser<GuesserType, typename StepperType::independent_variable_type, StateType>::value
  && std::is_same<IndVarType, typename StepperType::independent_variable_type>::value
  && std::is_same<StateType, typename StepperType::state_type>::value
  >
advance_n_steps_with_pre_step_guesser(StepperType & stepper,
				      StateType & state,
				      const IndVarType & start_val,
				      const IndVarType & step_size,
				      StepCount num_steps,
				      GuesserType && guesser)
{

  using observer_t = impl::NoOpStateObserver<IndVarType, StateType>;
  impl::advance_n_steps_with_fixed_dt(stepper, num_steps, start_val,
				      step_size, state,
				      observer_t(),
				      std::forward<GuesserType>(guesser));
}

template<
  class StepperType, class StateType, class StepSizePolicyType,
  class GuesserType, class IndVarType
  >
mpl::enable_if_t<
  Steppable<StepperType>::value
  && StepSizePolicy<StepSizePolicyType, typename StepperType::independent_variable_type>::value
  && StateGuesser<GuesserType, typename StepperType::independent_variable_type, StateType>::value
  && std::is_same<IndVarType, typename StepperType::independent_variable_type>::value
  && std::is_same<StateType, typename StepperType::state_type>::value
  >
advance_n_steps_with_pre_step_guesser(StepperType & stepper,
				      StateType & state,
				      const IndVarType & start_val,
				      StepSizePolicyType && step_size_policy,
				      StepCount num_steps,
				      GuesserType && guesser)
{

  using observer_t = impl::NoOpStateObserver<IndVarType, StateType>;
  impl::advance_n_steps_with_dt_policy(stepper, num_steps, start_val, state,
				       std::forward<StepSizePolicyType>(step_size_policy),
				       observer_t(),
				       std::forward<GuesserType>(guesser));
}

template<
  class StepperType, class StateType, class GuesserType,
  class ObserverType, class IndVarType
  >
mpl::enable_if_t<
  Steppable<StepperType>::value
  && StateGuesser<GuesserType, typename StepperType::independent_variable_type, StateType>::value
  && StateObserver<ObserverType, typename StepperType::independent_variable_type, StateType>::value
  && std::is_same<IndVarType, typename StepperType::independent_variable_type>::value
  && std::is_same<StateType, typename StepperType::state_type>::value
  >
advance_n_steps_with_pre_step_guesser(StepperType & stepper,
				      StateType & state,
				      const IndVarType & start_val,
				      const IndVarType & step_size,
				      StepCount num_steps,
				      GuesserType && guesser,
				      ObserverType && observer)
{

  impl::advance_n_steps_with_fixed_dt(stepper, num_steps, start_val,
				      step_size, state,
				      std::forward<ObserverType>(observer),
				      std::forward<GuesserType>(guesser));
}

template<
  class StepperType, class StateType, class StepSizePolicyType,
  class GuesserType, class ObserverType, class IndVarType
  >
mpl::enable_if_t<
  Steppable<StepperType>::value
  && StepSizePolicy<StepSizePolicyType, typename StepperType::independent_variable_type>::value
  && StateGuesser<GuesserType, typename StepperType::independent_variable_type, StateType>::value
  && StateObserver<ObserverType, typename StepperType::independent_variable_type, StateType>::value
  && std::is_same<IndVarType, typename StepperType::independent_variable_type>::value
  && std::is_same<StateType, typename StepperType::state_type>::value
  >
advance_n_steps_with_pre_step_guesser(StepperType & stepper,
				      StateType & state,
				      const IndVarType & start_val,
				      StepSizePolicyType && step_size_policy,
				      StepCount num_steps,
				      GuesserType && guesser,
				      ObserverType && observer)
{

  impl::advance_n_steps_with_dt_policy(stepper, num_steps, start_val, state,
				       std::forward<StepSizePolicyType>(step_size_policy),
				       std::forward<ObserverType>(observer),
				       std::forward<GuesserType>(guesser));
}

// ---------------------------------
// overload set for variadic stepper
// ---------------------------------

template<
  class StepperType, class StateType, class GuesserType, class IndVarType,
  class AuxT, class ...Args>
mpl::enable_if_t<
  SteppableWithAuxiliaryArgs<void, StepperType, AuxT, Args...>::value
  && StateGuesser<GuesserType, typename StepperType::independent_variable_type, StateType>::value
  && std::is_same<IndVarType, typename StepperType::independent_variable_type>::value
  && std::is_same<StateType, typename StepperType::state_type>::value
  >
advance_n_steps_with_pre_step_guesser(StepperType & stepper,
				      StateType & state,
				      const IndVarType & start_val,
				      const IndVarType & step_size,
				      StepCount num_steps,
				      GuesserType & guesser,
				      AuxT && auxArg,
				      Args && ... args)
{

  using observer_t = impl::NoOpStateObserver<IndVarType, StateType>;
  impl::advance_n_steps_with_fixed_dt(stepper, num_steps, start_val,
				      step_size, state,
				      observer_t(),
				      std::forward<GuesserType>(guesser),
				      std::forward<AuxT>(auxArg),
				      std::forward<Args>(args)...);
}

template<
  class StepperType, class StateType, class StepSizePolicyType,
  class GuesserType, class IndVarType, class AuxT, class ...Args>
mpl::enable_if_t<
  SteppableWithAuxiliaryArgs<void, StepperType, AuxT, Args...>::value
  && StepSizePolicy<StepSizePolicyType, typename StepperType::independent_variable_type>::value
  && StateGuesser<GuesserType, typename StepperType::independent_variable_type, StateType>::value
  && std::is_same<IndVarType, typename StepperType::independent_variable_type>::value
  && std::is_same<StateType, typename StepperType::state_type>::value
  >
advance_n_steps_with_pre_step_guesser(StepperType & stepper,
				      StateType & state,
				      const IndVarType & start_val,
				      StepSizePolicyType && step_size_policy,
				      StepCount num_steps,
				      GuesserType && guesser,
				      AuxT && auxArg,
				      Args && ... args)
{

  using observer_t = impl::NoOpStateObserver<IndVarType, StateType>;
  impl::advance_n_steps_with_dt_policy(stepper, num_steps, start_val, state,
				       std::forward<StepSizePolicyType>(step_size_policy),
				       observer_t(),
				       std::forward<GuesserType>(guesser),
				       std::forward<AuxT>(auxArg),
				       std::forward<Args>(args)...);
}

template<
  class StepperType, class StateType, class GuesserType,
  class ObserverType, class IndVarType, class AuxT, class ...Args
  >
mpl::enable_if_t<
  SteppableWithAuxiliaryArgs<void, StepperType, AuxT, Args...>::value
  && StateGuesser<GuesserType, typename StepperType::independent_variable_type, StateType>::value
  && StateObserver<ObserverType, typename StepperType::independent_variable_type, StateType>::value
  && std::is_same<IndVarType, typename StepperType::independent_variable_type>::value
  && std::is_same<StateType, typename StepperType::state_type>::value
  >
advance_n_steps_with_pre_step_guesser(StepperType & stepper,
				      StateType & state,
				      const IndVarType & start_val,
				      const IndVarType & step_size,
				      StepCount num_steps,
				      GuesserType & guesser,
				      ObserverType && observer,
				      AuxT && auxArg,
				      Args && ... args)
{

  impl::advance_n_steps_with_fixed_dt(stepper, num_steps, start_val,
				      step_size, state,
				      std::forward<ObserverType>(observer),
				      std::forward<GuesserType>(guesser),
				      std::forward<AuxT>(auxArg),
				      std::forward<Args>(args)...);
}

template<
  class StepperType, class StateType, class StepSizePolicyType,
  class GuesserType, class ObserverType, class IndVarType,
  class AuxT, class ...Args>
mpl::enable_if_t<
  SteppableWithAuxiliaryArgs<void, StepperType, AuxT, Args...>::value
  && StepSizePolicy<StepSizePolicyType, typename StepperType::independent_variable_type>::value
  && StateGuesser<GuesserType, typename StepperType::independent_variable_type, StateType>::value
  && StateObserver<ObserverType, typename StepperType::independent_variable_type, StateType>::value
  && std::is_same<IndVarType, typename StepperType::independent_variable_type>::value
  && std::is_same<StateType, typename StepperType::state_type>::value
  >
advance_n_steps_with_pre_step_guesser(StepperType & stepper,
				      StateType & state,
				      const IndVarType & start_val,
				      StepSizePolicyType && step_size_policy,
				      StepCount num_steps,
				      GuesserType && guesser,
				      ObserverType && observer,
				      AuxT && auxArg,
				      Args && ... args)
{

  impl::advance_n_steps_with_dt_policy(stepper, num_steps, start_val, state,
				       std::forward<StepSizePolicyType>(step_size_policy),
				       std::forward<ObserverType>(observer),
				       std::forward<GuesserType>(guesser),
				       std::forward<AuxT>(auxArg),
				       std::forward<Args>(args)...);
}

}}//end namespace pressio::ode
#endif