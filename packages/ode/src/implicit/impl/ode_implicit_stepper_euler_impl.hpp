/*
//@HEADER
// ************************************************************************
//
// ode_implicit_stepper_euler_impl.hpp
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

#ifndef ODE_IMPLICIT_IMPL_ODE_IMPLICIT_STEPPER_EULER_IMPL_HPP_
#define ODE_IMPLICIT_IMPL_ODE_IMPLICIT_STEPPER_EULER_IMPL_HPP_

namespace pressio{ namespace ode{ namespace implicitmethods{ namespace impl{

template<
  typename scalar_t,
  typename state_t,
  typename residual_t,
  typename jacobian_t,
  typename system_type,
  typename residual_policy_t,
  typename jacobian_policy_t
  >
class StepperBDF1
{
public:
  // these aliases are detected by solver
  using scalar_type	= scalar_t;
  using state_type	= state_t;
  using residual_type	= residual_t;
  using jacobian_type	= jacobian_t;

  static constexpr bool is_implicit = true;
  static constexpr bool is_explicit = false;
  static constexpr types::stepper_order_t order_value = 1;
  static constexpr std::size_t numAuxStates = 1;
  using tag_name = ::pressio::ode::implicitmethods::Euler;

  using aux_states_t =
    ::pressio::ode::AuxStatesManager<state_t, numAuxStates>;
  using standard_res_policy_t =
    ::pressio::ode::implicitmethods::policy::ResidualStandardPolicy<
    state_t, residual_t>;
  using standard_jac_policy_t =
    ::pressio::ode::implicitmethods::policy::JacobianStandardPolicy<
    state_t, jacobian_t>;

public:
  StepperBDF1() = delete;
  StepperBDF1(const StepperBDF1 & other)  = default;
  StepperBDF1 & operator=(const StepperBDF1 & other) = delete;
  StepperBDF1(StepperBDF1 && other)  = default;
  StepperBDF1 & operator=(StepperBDF1 && other) = delete;
  ~StepperBDF1() = default;

  StepperBDF1(const state_type & state,
	      const system_type & systemObj,
	      const residual_policy_t & resPolicyObj,
	      const jacobian_policy_t  & jacPolicyObj)
    : systemObj_{systemObj},
      auxStates_{state},
      residual_obj_{resPolicyObj},
      jacobian_obj_{jacPolicyObj}
  {}

  // cstr for standard residual and jacob policies
  template <
    typename T1 = residual_policy_t,
    typename T2 = jacobian_policy_t,
    ::pressio::mpl::enable_if_t<
      mpl::is_same<T1, standard_res_policy_t>::value and
      mpl::is_same<T2, standard_jac_policy_t>::value,
      int> = 0
    >
  StepperBDF1(const state_type & state,
	      const system_type & systemObj)
    : systemObj_{systemObj},
      auxStates_{state}
  {}

public:
  ::pressio::ode::types::stepper_order_t order() const
  {
    return order_value;
  }

  template<typename solver_type>
  void doStep(state_type & odeState,
	      const scalar_t & time,
	      const scalar_t & dt,
	      const types::step_t & step,
	      solver_type & solver)
  {
    PRESSIOLOG_DEBUG("bdf1 stepper: do step");

    static_assert(::pressio::ode::concepts::legitimate_solver_for_implicit_stepper<
      solver_type, decltype(*this), state_type>::value,
      "Invalid solver for BDF1 stepper");

    using nm1 = ode::nMinusOne;
    auto & odeState_nm1 = this->auxStates_.stateAt(nm1());
    this->dt_ = dt;
    this->t_ = time;
    this->step_ = step;
    ::pressio::ops::deep_copy(odeState_nm1, odeState);

    try{
      solver.solve(*this, odeState);
    }
    catch (::pressio::eh::nonlinear_solve_failure const & e)
    {
      // the state before attempting solution was stored in y_n-1,
      // so revert odeState to that
      auto & rollBackState = this->auxStates_.stateAt(ode::nMinusOne());
      ::pressio::ops::deep_copy(odeState, rollBackState);

      // now throw
      throw ::pressio::eh::time_step_failure();
    }
  }

  template<typename solver_type, typename guess_callback_t>
  void doStep(state_type & odeState,
	      const scalar_t & time,
	      const scalar_t & dt,
	      const types::step_t & step,
	      solver_type & solver,
	      guess_callback_t && guesserCb)
  {
    PRESSIOLOG_DEBUG("bdf1 stepper: do step with callback to state guesser");

    static_assert(::pressio::ode::concepts::legitimate_solver_for_implicit_stepper<
      solver_type, decltype(*this), state_type>::value,
      "Invalid solver for BDF1 stepper");

    using nm1 = ode::nMinusOne;
    auto & odeState_nm1 = this->auxStates_.stateAt(nm1());
    this->dt_ = dt;
    this->t_ = time;
    this->step_ = step;
    ::pressio::ops::deep_copy(odeState_nm1, odeState);
    guesserCb(step, time, odeState);


    try{
      solver.solve(*this, odeState);
    }
    catch (::pressio::eh::nonlinear_solve_failure const & e)
    {
      // the state before attempting solution was stored in y_n-1,
      // so revert odeState to that
      auto & rollBackState = this->auxStates_.stateAt(ode::nMinusOne());
      ::pressio::ops::deep_copy(odeState, rollBackState);

      // now throw
      throw ::pressio::eh::time_step_failure();
    }
  }

  residual_t createResidual() const
  {
    const auto & resPol = static_cast<const residual_policy_t&>(residual_obj_);
    return resPol.create(systemObj_.get());
  }

  jacobian_t createJacobian() const
  {
    const auto & jacPol = static_cast<const jacobian_policy_t&>(jacobian_obj_);
    return jacPol.create(systemObj_.get());
  }

  void residual(const state_t & odeState, residual_t & R) const
  {
    const auto & resPol = static_cast<const residual_policy_t&>(residual_obj_);
    resPol.template compute<tag_name>(
      odeState, this->auxStates_, this->systemObj_.get(),
      this->t_, this->dt_, this->step_, R);
  }

  void jacobian(const state_t & odeState, jacobian_t & J) const
  {
    const auto & jacPol = static_cast<const jacobian_policy_t&>(jacobian_obj_);
    jacPol.template compute<
      tag_name>(odeState, this->auxStates_, this->systemObj_.get(),
                this->t_, this->dt_, this->step_, J);
  }

private:
  scalar_t t_  = {};
  scalar_t dt_ = {};
  types::step_t step_  = {};
  std::reference_wrapper<const system_type> systemObj_;
  aux_states_t auxStates_;

  // conditional type of the policies
  // if we have a standard policy, then it is default constr
  // if we have a user-defined policy, we take a const & to it
  typename std::conditional<
    mpl::is_same<standard_res_policy_t, residual_policy_t>::value,
    const residual_policy_t,
    std::reference_wrapper<const residual_policy_t>
    >::type residual_obj_;

  typename std::conditional<
    mpl::is_same<standard_jac_policy_t, jacobian_policy_t>::value,
    const jacobian_policy_t,
    std::reference_wrapper<const jacobian_policy_t>
    >::type jacobian_obj_;
};

}}}}
#endif  // ODE_IMPLICIT_IMPL_ODE_IMPLICIT_STEPPER_EULER_IMPL_HPP_
