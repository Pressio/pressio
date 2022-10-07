/*
//@HEADER
// ************************************************************************
//
// ode_explicit_compose.hpp
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

#ifndef ODE_STEPPERS_IMPL_ODE_EXPLICIT_CREATE_STEPPER_IMPL_HPP_
#define ODE_STEPPERS_IMPL_ODE_EXPLICIT_CREATE_STEPPER_IMPL_HPP_

#include "ode_explicit_stepper_without_mass_matrix.hpp"
#include "ode_explicit_stepper_with_mass_matrix.hpp"

namespace pressio{ namespace ode{ namespace impl{

template<class ImplClassType, class SystemType>
auto create_explicit_stepper(StepScheme name,
			     SystemType && system)
{

  if (name == StepScheme::ForwardEuler){
    return ImplClassType(ode::ForwardEuler(),
			 std::forward<SystemType>(system));
  }

  else if (name == StepScheme::RungeKutta4){
    return ImplClassType(ode::RungeKutta4(),
			 std::forward<SystemType>(system));
  }

  else if (name == StepScheme::AdamsBashforth2){
    return ImplClassType(ode::AdamsBashforth2(),
			 std::forward<SystemType>(system));
  }

  else if (name == StepScheme::SSPRungeKutta3){
    return ImplClassType(ode::SSPRungeKutta3(),
			 std::forward<SystemType>(system));
  }

  else{
    throw std::runtime_error("ode:: create_explicit_stepper: invalid StepScheme enum value");
  }

}

}}}
#endif  // ODE_STEPPERS_IMPL_ODE_EXPLICIT_COMPOSE_HPP_
