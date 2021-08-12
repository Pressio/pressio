/*
//@HEADER
// ************************************************************************
//
// ode_implicit_stepper.hpp
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

#ifndef ODE_IMPLICIT_ODE_IMPLICIT_STEPPER_HPP_
#define ODE_IMPLICIT_ODE_IMPLICIT_STEPPER_HPP_

#include "./impl/ode_implicit_stepper_compose.hpp"

namespace pressio{ namespace ode{

// BDF1
template<class SystemType, class StateType, class ...Args>
auto create_bdf1_stepper(const SystemType & system,  const StateType & state,  Args && ... args)
->  decltype( 
         impl::create_stepper_impl<implicitmethods::BDF1>(system, state, std::forward<Args>(args)...) )
{
  return impl::create_stepper_impl<implicitmethods::BDF1>(system, state, std::forward<Args>(args)...);
};

template<class ResidualType, class JacobianType, class SystemType, class StateType, class ...Args>
auto create_bdf1_stepper_partial_deduction(const SystemType & system, const StateType & state, Args && ... args)
->  decltype( 
         impl::create_stepper_partial_deduction_impl<
  implicitmethods::BDF1, ResidualType, JacobianType>(system, state, std::forward<Args>(args)...) )

{
  return impl::create_stepper_partial_deduction_impl<
  implicitmethods::BDF1, ResidualType, JacobianType>(system, state, std::forward<Args>(args)...);
};

// BDF2
template<class SystemType, class StateType, class ...Args>
auto create_bdf2_stepper(const SystemType & system,  const StateType & state,  Args && ... args)
->  decltype( 
         impl::create_stepper_impl<implicitmethods::BDF2>(system, state, std::forward<Args>(args)...) )
{
  return impl::create_stepper_impl<implicitmethods::BDF2>(system, state, std::forward<Args>(args)...);
};

template<class ResidualType, class JacobianType, class SystemType, class StateType, class ...Args>
auto create_bdf2_stepper_partial_deduction(const SystemType & system, const StateType & state, Args && ... args)
->  decltype( 
         impl::create_stepper_partial_deduction_impl<
  implicitmethods::BDF2, ResidualType, JacobianType>(system, state, std::forward<Args>(args)...) )

{
  return impl::create_stepper_partial_deduction_impl<
  implicitmethods::BDF2, ResidualType, JacobianType>(system, state, std::forward<Args>(args)...);
};

// CrankNicolson
template<class SystemType, class StateType, class ...Args>
auto create_cranknicolson_stepper(const SystemType & system,  const StateType & state,  Args && ... args)
->  decltype(   
         impl::create_stepper_impl<implicitmethods::CrankNicolson>(system, state, std::forward<Args>(args)...) )
{
  return impl::create_stepper_impl<implicitmethods::CrankNicolson>(system, state, std::forward<Args>(args)...);
};

template<class ResidualType, class JacobianType, class SystemType, class StateType, class ...Args>
auto create_cranknicolson_stepper_partial_deduction(const SystemType & system, const StateType & state, Args && ... args)
->  decltype( 
         impl::create_stepper_partial_deduction_impl<
  implicitmethods::CrankNicolson, ResidualType, JacobianType>(system, state, std::forward<Args>(args)...) )

{
  return impl::create_stepper_partial_deduction_impl<
  implicitmethods::CrankNicolson, ResidualType, JacobianType>(system, state, std::forward<Args>(args)...);
};


//
// Arbitrary
//
template<
  int order, 
  int num_states, 
  class SystemType, 
  class StateType, 
  class ...Args,
  class ReturnType = impl::ImplicitCompose_t<
    implicitmethods::Arbitrary, 
    StepperOrder<order>, StepperTotalNumberOfStates<num_states>, SystemType, StateType, Args...>
  >
ReturnType create_arbitrary_stepper(const SystemType & system, const StateType & state, Args && ... args)
{
  return ReturnType(state, system, std::forward<Args>(args)...);
};

template<
  int order, 
  int num_states, 
  class ResidualType, 
  class JacobianType, 
  class SystemType, 
  class StateType, 
  class ReturnType = impl::ImplicitCompose_t<
    implicitmethods::Arbitrary,
    StepperOrder<order>, StepperTotalNumberOfStates<num_states>,
    SystemType, StateType, ResidualType, JacobianType, void>
  >
ReturnType create_arbitrary_stepper_partial_deduction(const SystemType & system, const StateType & state)
{
  return ReturnType(state, system);
};


template<
  int order,
  int num_states,
  class ResidualType,
  class JacobianType,
  class SystemType,
  class StateType,
  class ResidualPolicyType,
  class JacobianPolicyType,
  class ReturnType = impl::ImplicitCompose_t<
    implicitmethods::Arbitrary,
    StepperOrder<order>, StepperTotalNumberOfStates<num_states>,
    SystemType, StateType, ResidualType, JacobianType, ResidualPolicyType, JacobianPolicyType>
  >
ReturnType create_arbitrary_stepper_partial_deduction(const SystemType & system, const StateType & state, 
                                                      ResidualPolicyType && rPol, JacobianPolicyType && jPol)
{
  return ReturnType(state, system,
		    std::forward<ResidualPolicyType>(rPol),
		    std::forward<JacobianPolicyType>(jPol));
};

}} // end namespace pressio::ode
#endif  // ODE_IMPLICIT_ODE_IMPLICIT_STEPPER_HPP_
