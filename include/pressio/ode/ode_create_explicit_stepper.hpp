/*
//@HEADER
// ************************************************************************
//
// ode_create_explicit_stepper.hpp
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

#ifndef PRESSIO_ODE_ODE_CREATE_EXPLICIT_STEPPER_HPP_
#define PRESSIO_ODE_ODE_CREATE_EXPLICIT_STEPPER_HPP_

#include "./impl/ode_explicit_stepper_without_mass_matrix.hpp"
#include "./impl/ode_explicit_stepper_with_mass_matrix.hpp"
#include "./impl/ode_explicit_create_impl.hpp"

namespace pressio{ namespace ode{

//
// basic, no mass matrix
//
template<
  class SystemType,
  std::enable_if_t<
    RealValuedOdeSystem<mpl::remove_cvref_t<SystemType>>::value,
    int > = 0
  >
auto create_explicit_stepper(StepScheme schemeName,                     // (1)
			     SystemType && odeSystem)
{

  using sys_type = mpl::remove_cvref_t<SystemType>;
  using ind_var_type = typename sys_type::independent_variable_type;
  using state_type   = typename sys_type::state_type;
  using rhs_type = typename sys_type::rhs_type;

  /* IMPORTANT: use "SystemType" as template arg because that it the
     right type carrying how we store the system and NOT SystemType &&
     for the following reason: when user passes a non-temporary object,
     SystemType is deduced to be a reference, so the concrete stepper class
     will hold a reference to the provided system object.
     When the user passes system to be a temporary object,
     SystemType will be deduced so that the stepper will hold an **instance**
     of the system that is move-constructed (if applicable, or copy-constructed)
     from the system argument.
  */
  using impl_type = impl::ExplicitStepperNoMassMatrixImpl<
    state_type, ind_var_type, SystemType, rhs_type>;
  return impl::create_explicit_stepper<impl_type>
    (schemeName, std::forward<SystemType>(odeSystem));
}

//
// WITH mass matrix
//
template<
  class SystemType,
  std::enable_if_t<
    RealValuedOdeSystemFusingMassMatrixAndRhs<mpl::remove_cvref_t<SystemType>>::value,
    int > = 0
  >
auto create_explicit_stepper(StepScheme schemeName,                     // (2)
			     SystemType && odeSystem)
{

  using sys_type = mpl::remove_cvref_t<SystemType>;
  using ind_var_type = typename sys_type::independent_variable_type;
  using state_type   = typename sys_type::state_type;
  using rhs_type = typename sys_type::rhs_type;

  // use "SystemType" as template arg, see above for reason
  using impl_type = impl::ExplicitStepperWithMassMatrixImpl<
    state_type, ind_var_type, SystemType, rhs_type>;
  return impl::create_explicit_stepper<impl_type>
    (schemeName, std::forward<SystemType>(odeSystem));
}

//
// auxiliary scheme-specific functions
//
template<class ...Args>
auto create_forward_euler_stepper(Args && ...args){
  return create_explicit_stepper(StepScheme::ForwardEuler,
				 std::forward<Args>(args)...);
}

template<class ...Args>
auto create_rk4_stepper(Args && ...args){
  return create_explicit_stepper(StepScheme::RungeKutta4,
				 std::forward<Args>(args)...);
}

template<class ...Args>
auto create_ab2_stepper(Args && ...args){
  return create_explicit_stepper(StepScheme::AdamsBashforth2,
				 std::forward<Args>(args)...);
}

template<class ...Args>
auto create_ssprk3_stepper(Args && ...args){
  return create_explicit_stepper
    (StepScheme::SSPRungeKutta3, std::forward<Args>(args)...);
}

}} // end namespace pressio::ode
#endif  // PRESSIO_ODE_ODE_CREATE_EXPLICIT_STEPPER_HPP_
