/*
//@HEADER
// ************************************************************************
//
// ode_explicit_stepper_base.hpp
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

#ifndef ODE_STEPPERS_EXPLICIT_STEPPERS_BASE_EXPLICIT_STEPPER_BASE_HPP_
#define ODE_STEPPERS_EXPLICIT_STEPPERS_BASE_EXPLICIT_STEPPER_BASE_HPP_

namespace pressio{ namespace ode{ namespace explicitmethods{

/*
 * (1) constructors here should be private but we need
 * them public to enable interfacing with pybind11
 */

template<typename stepper_type>
class StepperBase
{
private:
  using stepper_traits	= ::pressio::ode::details::traits<stepper_type>;
  using scalar_t	= typename stepper_traits::scalar_t;
  using state_t		= typename stepper_traits::state_t;
  // friend the derived so that it can access private constructors
  friend stepper_type;

public:
  typename stepper_traits::order_t order() const{
    return stepper_traits::order_value;
  }

  void operator()(state_t & odeStateInOut,
  		  const scalar_t & time,
  		  const scalar_t & dt,
  		  const ::pressio::ode::types::step_t & step)
  {
    static_cast<stepper_type&>(*this).doStep(odeStateInOut, time, dt, step);
  }

private:
  StepperBase()  = default;
  ~StepperBase() = default;

  // copy cnstr
  StepperBase(const StepperBase & other)  = delete;
  // copy assignment
  StepperBase & operator=(const StepperBase & other)  = delete;
  // move cnstr
  StepperBase(StepperBase && other)  = delete;
  // move assign
  StepperBase & operator=(StepperBase && other)  = delete;

};//end class

}}}//end namespace pressio::ode::explicitmethods
#endif
