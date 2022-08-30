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

#ifndef ODE_ADVANCERS_ODE_ADVANCE_N_STEPS_HPP_
#define ODE_ADVANCERS_ODE_ADVANCE_N_STEPS_HPP_

#include "./impl/ode_advance_noop_observer.hpp"
#include "./impl/ode_advance_noop_guesser.hpp"
#include "./impl/ode_advance_n_steps.hpp"
#include "./impl/ode_mandates.hpp"

namespace pressio{ namespace ode{

// ---------------------------------
// overload set for steppable
// ---------------------------------

// const dt
template<
  class StepperType, class StateType, class IndVarType
  >
mpl::enable_if_t<
  Steppable<StepperType>::value
  >
advance_n_steps(StepperType & stepper,
		StateType & state,
		const IndVarType & start_val,
		const IndVarType & step_size,
		StepCount num_steps)
{

  impl::mandate_on_ind_var_and_state_types(stepper, state, start_val);
  using observer_t = impl::NoOpStateObserver<IndVarType, StateType>;
  using guesser_t  = impl::NoOpStateGuesser<IndVarType, StateType>;
  impl::advance_n_steps_with_fixed_dt(stepper, num_steps, start_val,
				      step_size, state,
				      observer_t(), guesser_t());
}

// dt policy provided
template<
  class StepperType, class StateType,
  class StepSizePolicyType, class IndVarType
  >
mpl::enable_if_t<
  Steppable<StepperType>::value
  && StepSizePolicy<StepSizePolicyType, typename StepperType::independent_variable_type>::value
  >
advance_n_steps(StepperType & stepper,
		StateType & state,
		const IndVarType & start_val,
		StepSizePolicyType && step_size_policy,
		StepCount num_steps)
{

  impl::mandate_on_ind_var_and_state_types(stepper, state, start_val);
  using observer_t = impl::NoOpStateObserver<IndVarType, StateType>;
  using guesser_t  = impl::NoOpStateGuesser<IndVarType, StateType>;
  impl::advance_n_steps_with_dt_policy(stepper, num_steps, start_val, state,
				       std::forward<StepSizePolicyType>(step_size_policy),
				       observer_t(), guesser_t());
}

// const dt and observer
template<
  class StepperType, class StateType, class ObserverType, class IndVarType
  >
mpl::enable_if_t<
  Steppable<StepperType>::value
  && StateObserver<ObserverType, typename StepperType::independent_variable_type, StateType>::value
  >
advance_n_steps(StepperType & stepper,
		StateType & state,
		const IndVarType & start_val,
		const IndVarType & step_size,
		StepCount num_steps,
		ObserverType && observer)
{

  impl::mandate_on_ind_var_and_state_types(stepper, state, start_val);
  using guesser_t  = impl::NoOpStateGuesser<IndVarType, StateType>;
  impl::advance_n_steps_with_fixed_dt(stepper, num_steps, start_val,
				      step_size, state,
				      std::forward<ObserverType>(observer),
				      guesser_t());
}

// dt policy provided and observer
template<
  class StepperType, class StateType, class ObserverType,
  class StepSizePolicyType, class IndVarType
  >
mpl::enable_if_t<
  Steppable<StepperType>::value
  && StateObserver<ObserverType, typename StepperType::independent_variable_type, StateType>::value
  && StepSizePolicy<StepSizePolicyType, typename StepperType::independent_variable_type>::value
  >
advance_n_steps(StepperType & stepper,
		StateType & state,
		const IndVarType & start_val,
		StepSizePolicyType && step_size_policy,
		StepCount num_steps,
		ObserverType && observer)
{

  impl::mandate_on_ind_var_and_state_types(stepper, state, start_val);
  using guesser_t  = impl::NoOpStateGuesser<IndVarType, StateType>;
  impl::advance_n_steps_with_dt_policy(stepper, num_steps, start_val, state,
				       std::forward<StepSizePolicyType>(step_size_policy),
				       std::forward<ObserverType>(observer),
				       guesser_t());
}


// ---------------------------------
// overload set for steppable with extra args
// ---------------------------------

// const dt
template<
  class StepperType, class StateType, class IndVarType,
  class AuxT, class ...Args>
mpl::enable_if_t<
  SteppableWithAuxiliaryArgs<void, StepperType, AuxT, Args...>::value
  && !StateObserver<AuxT, typename StepperType::independent_variable_type, StateType>::value
  >
advance_n_steps(StepperType & stepper,
		StateType & state,
		const IndVarType & start_val,
		const IndVarType & step_size,
		StepCount num_steps,
		AuxT && auxArg,
		Args && ... args)
{

  impl::mandate_on_ind_var_and_state_types(stepper, state, start_val);
  using observer_t = impl::NoOpStateObserver<IndVarType, StateType>;
  using guesser_t  = impl::NoOpStateGuesser<IndVarType, StateType>;
  observer_t observer;
  impl::advance_n_steps_with_fixed_dt(stepper, num_steps, start_val,
				      step_size, state, observer,
				      guesser_t(),
				      std::forward<AuxT>(auxArg),
				      std::forward<Args>(args)...);
}

// dt policy provided
template<
  class StepperType, class StateType, class StepSizePolicyType, class IndVarType,
  class AuxT, class ...Args>
mpl::enable_if_t<
  SteppableWithAuxiliaryArgs<void, StepperType, AuxT, Args...>::value
  && StepSizePolicy<StepSizePolicyType, typename StepperType::independent_variable_type>::value
  && !StateObserver<AuxT, typename StepperType::independent_variable_type, StateType>::value
  >
advance_n_steps(StepperType & stepper,
		StateType & state,
		const IndVarType & start_val,
		StepSizePolicyType && step_size_policy,
		StepCount num_steps,
		AuxT && auxArg,
		Args && ... args)
{

  impl::mandate_on_ind_var_and_state_types(stepper, state, start_val);
  using observer_t = impl::NoOpStateObserver<IndVarType, StateType>;
  using guesser_t  = impl::NoOpStateGuesser<IndVarType, StateType>;
  impl::advance_n_steps_with_dt_policy(stepper, num_steps, start_val, state,
				       std::forward<StepSizePolicyType>(step_size_policy),
				       observer_t(), guesser_t(),
				       std::forward<AuxT>(auxArg),
				       std::forward<Args>(args)...);
}

// const dt and observer
template<
  class StepperType, class StateType, class ObserverType, class IndVarType,
  class AuxT, class ...Args>
mpl::enable_if_t<
  SteppableWithAuxiliaryArgs<void, StepperType, AuxT, Args...>::value
  && StateObserver<ObserverType, typename StepperType::independent_variable_type, StateType>::value
  && !StateObserver<AuxT, typename StepperType::independent_variable_type, StateType>::value
  >
advance_n_steps(StepperType & stepper,
		StateType & state,
		const IndVarType & start_val,
		const IndVarType & step_size,
		StepCount num_steps,
		ObserverType && observer,
		AuxT && auxArg,
		Args && ... args)
{

  impl::mandate_on_ind_var_and_state_types(stepper, state, start_val);
  using guesser_t  = impl::NoOpStateGuesser<IndVarType, StateType>;
  impl::advance_n_steps_with_fixed_dt(stepper, num_steps, start_val,
				      step_size, state,
				      std::forward<ObserverType>(observer),
				      guesser_t(),
				      std::forward<AuxT>(auxArg),
				      std::forward<Args>(args)...);
}

// dt policy provided and observer
template<
  class StepperType, class StateType, class StepSizePolicyType,
  class ObserverType, class IndVarType,
  class AuxT, class ...Args>
mpl::enable_if_t<
  SteppableWithAuxiliaryArgs<void, StepperType, AuxT, Args...>::value
  && StepSizePolicy<StepSizePolicyType, typename StepperType::independent_variable_type>::value
  && StateObserver<ObserverType, typename StepperType::independent_variable_type, StateType>::value
  && !StateObserver<AuxT, typename StepperType::independent_variable_type, StateType>::value
  >
advance_n_steps(StepperType & stepper,
		StateType & state,
		const IndVarType & start_val,
		StepSizePolicyType && step_size_policy,
		StepCount num_steps,
		ObserverType && observer,
		AuxT && auxArg,
		Args && ... args)
{

  impl::mandate_on_ind_var_and_state_types(stepper, state, start_val);
  using guesser_t  = impl::NoOpStateGuesser<IndVarType, StateType>;
  impl::advance_n_steps_with_dt_policy(stepper, num_steps, start_val, state,
				       std::forward<StepSizePolicyType>(step_size_policy),
				       std::forward<ObserverType>(observer),
				       guesser_t(),
				       std::forward<AuxT>(auxArg),
				       std::forward<Args>(args)...);
}

}} //end namespace pressio::ode
#endif
