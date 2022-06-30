/*
//@HEADER
// ************************************************************************
//
// ode_implicit_policy_jacobian.hpp
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

#ifndef ODE_STEPPERS_IMPL_ODE_IMPLICIT_POLICY_JACOBIAN_HPP_
#define ODE_STEPPERS_IMPL_ODE_IMPLICIT_POLICY_JACOBIAN_HPP_

namespace pressio{ namespace ode{ namespace impl{

template<
  class SystemType,
  class IndVarType,
  class StateType,
  class JacobianType
  >
class JacobianStandardPolicy
{
public:
  // required
  using independent_variable_type = IndVarType;
  using state_type    = StateType;
  using jacobian_type = JacobianType;

public:
  JacobianStandardPolicy() = delete;

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  explicit JacobianStandardPolicy(SystemType systemIn)
    : systemObj_(systemIn){}
#else
  explicit JacobianStandardPolicy(SystemType && systemIn)
    : systemObj_( std::forward<SystemType>(systemIn) ){}
#endif

  JacobianStandardPolicy(const JacobianStandardPolicy &) = default;
  JacobianStandardPolicy & operator=(const JacobianStandardPolicy &) = default;
  JacobianStandardPolicy(JacobianStandardPolicy &&) = default;
  JacobianStandardPolicy & operator=(JacobianStandardPolicy &&) = default;
  ~JacobianStandardPolicy() = default;

public:
  StateType createState() const{
    StateType result(systemObj_.get().createState());
    return result;
  }

  JacobianType createJacobian() const{
    JacobianType JJ(systemObj_.get().createJacobian());
    return JJ;
  }

  template <class StencilStatesContainerType>
  void operator()(StepScheme name,
		  const StateType & odeCurrentState,
		  const StencilStatesContainerType & stencilStates,
		  const ::pressio::ode::StepEndAt<IndVarType> & evalTime,
		  ::pressio::ode::StepCount step,
		  const ::pressio::ode::StepSize<IndVarType> & dt,
		  JacobianType & J) const
  {
    systemObj_.get().jacobian(odeCurrentState, evalTime.get(), J);

    if (name == StepScheme::BDF1){
      ::pressio::ode::impl::discrete_jacobian(ode::BDF1(), J, dt.get());
    }
    else if (name == StepScheme::BDF2){
      ::pressio::ode::impl::discrete_jacobian(ode::BDF2(), J, dt.get());
    }
    else if (name == StepScheme::CrankNicolson){
      ::pressio::ode::impl::discrete_jacobian(ode::CrankNicolson(), J, dt.get());
    }
  }

private:
  ::pressio::utils::InstanceOrReferenceWrapper<SystemType> systemObj_;
};

}}}//end namespace pressio::ode::implicitmethods::policy
#endif  // ODE_STEPPERS_IMPL_ODE_IMPLICIT_POLICY_JACOBIAN_HPP_
