/*
//@HEADER
// ************************************************************************
//
// ode_implicit_stepper_arbitrary_impl.hpp
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

#ifndef ODE_IMPLICIT_IMPL_ODE_IMPLICIT_STEPPER_ARBITRARY_IMPL_HPP_
#define ODE_IMPLICIT_IMPL_ODE_IMPLICIT_STEPPER_ARBITRARY_IMPL_HPP_

namespace pressio{ namespace ode{ namespace implicitmethods{ namespace impl{

template<
  typename scalar_t,
  typename ode_state_type,
  typename ode_residual_type,
  typename ode_jacobian_type,
  typename system_type,
  typename order_setter_t,
  typename tot_n_setter_t,
  typename residual_policy_t,
  typename jacobian_policy_t
  >
class StepperArbitrary
{
public:
  // these need to be here because are detected by solver
  using scalar_type	= scalar_t;
  using state_type	= ode_state_type;
  using residual_type	= ode_residual_type;
  using jacobian_type	= ode_jacobian_type;


  static constexpr bool is_implicit = true;
  static constexpr bool is_explicit = false;
  // numAuxStates is the number of auxiliary states needed, so all other beside y_n
  static constexpr std::size_t numAuxStates = tot_n_setter_t::value - 1;
  using tag_name = ::pressio::ode::implicitmethods::Arbitrary;

  using aux_states_t =
    ::pressio::ode::AuxStatesManager<ode_state_type, numAuxStates>;

  using standard_res_policy_t =
    ::pressio::ode::implicitmethods::policy::ResidualStandardPolicy<
    ode_state_type, ode_residual_type>;
  using standard_jac_policy_t =
    ::pressio::ode::implicitmethods::policy::JacobianStandardPolicy<
    ode_state_type, ode_jacobian_type>;

public:
  StepperArbitrary() = delete;
  StepperArbitrary(const StepperArbitrary & other)  = default;
  StepperArbitrary & operator=(const StepperArbitrary & other) = delete;
  StepperArbitrary(StepperArbitrary && other)  = default;
  StepperArbitrary & operator=(StepperArbitrary && other) = delete;
  ~StepperArbitrary() = default;

  // note that here residual_policy_t can be a reference already
  // so we don't need to specify & in argument to constructor
  StepperArbitrary(const ode_state_type & state,
		   const system_type & systemObj,
		   residual_policy_t resPolicyObj,
		   jacobian_policy_t jacPolicyObj)
    : systemObj_{systemObj},
      auxStates_{state},
      recoveryState_{state},
      resPolicy_{resPolicyObj},
      jacPolicy_{jacPolicyObj}
    {}

  // cstr for standard residual and jacob policies
  template <
    typename T1 = residual_policy_t,
    typename T2 = jacobian_policy_t,
    ::pressio::mpl::enable_if_t<
      std::is_default_constructible<T1>::value and
      std::is_default_constructible<T2>::value,
      int > = 0
    >
  StepperArbitrary(const ode_state_type & state,
                   const system_type & systemObj)
    : systemObj_{systemObj},
      auxStates_{state},
      recoveryState_{state},
      resPolicy_{T1()},
      jacPolicy_{T2()}
  {}

public:
  ::pressio::ode::types::stepper_order_t order() const
  {
    return order_setter_t::value;
  }

  residual_type createResidual() const
  {
    return resPolicy_.get().create(systemObj_.get());
  }

  jacobian_type createJacobian() const
  {
    return jacPolicy_.get().create(systemObj_.get());
  }

  void residual(const state_type & odeState, residual_type & R) const
  {
    resPolicy_.get().template compute
      <tag_name>(odeState, this->auxStates_,
		 this->systemObj_.get(),
		 this->t_, this->dt_, this->step_, R);
  }

  void jacobian(const state_type & odeState, jacobian_type & J) const
  {
    jacPolicy_.get().template compute
      <tag_name>(odeState, this->auxStates_,
		 this->systemObj_.get(),
		 this->t_, this->dt_, this->step_, J);
  }

  template<typename solver_type>
  void doStep(ode_state_type & odeState,
	      const scalar_t & t,
	      const scalar_t & dt,
	      const types::step_t & step,
	      solver_type & solver)
  {
    PRESSIOLOG_DEBUG("arbitrary stepper: do step");

    static_assert
      (::pressio::ode::constraints::legitimate_solver_for_implicit_stepper<
      solver_type, decltype(*this), ode_state_type>::value,
      "Invalid solver for Arbitrary stepper");

    this->dt_ = dt;
    this->t_ = t;
    this->step_ = step;

    constexpr auto nAux = decltype(this->auxStates_)::size();
    this->storeRecoveryState<nAux>();
    this->updateAuxiliaryStorage<nAux>(odeState);

    try{
      solver.solve(*this, odeState);
    }
    catch (::pressio::eh::nonlinear_solve_failure const & e)
    {
      // the state before attempting solution was stored in y_n-1,
      // so revert odeState to that
      auto & rollBackState = this->auxStates_.stateAt(ode::nMinusOne());
      ::pressio::ops::deep_copy(odeState, rollBackState);

      // we also need to revert back all the auxiliary states
      this->rollBackAuxiliaryStorage<nAux>();

      // now throw
      throw ::pressio::eh::time_step_failure();
    }
  }

private:
  //-------------------
  // one aux states
  //-------------------
  template<std::size_t nAux, mpl::enable_if_t<nAux==1, int > = 0>
  void storeRecoveryState()
  {
    auto & src = this->auxStates_.stateAt(ode::nMinusOne());
    ::pressio::ops::deep_copy(recoveryState_, src);
  }

  template<std::size_t nAux, mpl::enable_if_t<nAux==1, int > = 0>
  void updateAuxiliaryStorage(const ode_state_type & odeState)
  {
    // copy y_n into y_n-1
    auto & y_nm1 = this->auxStates_.stateAt(ode::nMinusOne());
    ::pressio::ops::deep_copy(y_nm1, odeState);
  }

  template<std::size_t nAux, mpl::enable_if_t<nAux==1, int > = 0>
  void rollBackAuxiliaryStorage()
  {
    auto & y_nm1 = this->auxStates_.stateAt(ode::nMinusOne());
    ::pressio::ops::deep_copy(y_nm1, recoveryState_);
  }

  //-------------------
  // two aux states
  //-------------------
  template<std::size_t nAux, mpl::enable_if_t<nAux==2, int > = 0>
  void storeRecoveryState()
  {
    auto & src = this->auxStates_.stateAt(ode::nMinusTwo());
    ::pressio::ops::deep_copy(recoveryState_, src);
  }

  template<std::size_t nAux, mpl::enable_if_t<nAux==2, int> = 0>
  void updateAuxiliaryStorage(const ode_state_type & odeState)
  {
    // copy y_n-1 into y_n-2
    auto & y_nm1 = this->auxStates_.stateAt(ode::nMinusOne());
    auto & y_nm2 = this->auxStates_.stateAt(ode::nMinusTwo());
    ::pressio::ops::deep_copy(y_nm2, y_nm1);
    // copy y_n into y_n-1
    ::pressio::ops::deep_copy(y_nm1, odeState);
  }

  template<std::size_t nAux, mpl::enable_if_t<nAux==2, int> = 0>
  void rollBackAuxiliaryStorage()
  {
    // copy y_n-2 into y_n-1
    auto & y_nm1 = this->auxStates_.stateAt(ode::nMinusOne());
    auto & y_nm2 = this->auxStates_.stateAt(ode::nMinusTwo());
    ::pressio::ops::deep_copy(y_nm1, y_nm2);
    //copy safestate to y_n-2
    ::pressio::ops::deep_copy(y_nm2, recoveryState_);
  }

  //-------------------
  // three aux states
  //-------------------
  template<std::size_t nAux, mpl::enable_if_t<nAux==3, int > = 0>
  void storeRecoveryState()
  {
    auto & src = this->auxStates_.stateAt(ode::nMinusThree());
    ::pressio::ops::deep_copy(recoveryState_, src);
  }

  template<std::size_t nAux, mpl::enable_if_t<nAux==3, int> = 0>
  void updateAuxiliaryStorage(const ode_state_type & odeState)
  {
    auto & y_nm1 = this->auxStates_.stateAt(ode::nMinusOne());
    auto & y_nm2 = this->auxStates_.stateAt(ode::nMinusTwo());
    auto & y_nm3 = this->auxStates_.stateAt(ode::nMinusThree());
    ::pressio::ops::deep_copy(y_nm3, y_nm2);
    ::pressio::ops::deep_copy(y_nm2, y_nm1);
    ::pressio::ops::deep_copy(y_nm1, odeState);
  }

  template<std::size_t nAux, mpl::enable_if_t<nAux==3, int> = 0>
  void rollBackAuxiliaryStorage()
  {
    auto & y_nm1 = this->auxStates_.stateAt(ode::nMinusOne());
    auto & y_nm2 = this->auxStates_.stateAt(ode::nMinusTwo());
    auto & y_nm3 = this->auxStates_.stateAt(ode::nMinusThree());
    ::pressio::ops::deep_copy(y_nm1, y_nm2);
    ::pressio::ops::deep_copy(y_nm2, y_nm3);
    ::pressio::ops::deep_copy(y_nm3, recoveryState_);
  }

private:
  scalar_t t_  = {};
  scalar_t dt_ = {};
  types::step_t step_  = {};
  std::reference_wrapper<const system_type> systemObj_;

  aux_states_t auxStates_;

  // state object to ensure the strong guarantee for handling excpetions
  ode_state_type recoveryState_;

  ::pressio::utils::instance_or_reference_wrapper<residual_policy_t> resPolicy_;
  ::pressio::utils::instance_or_reference_wrapper<jacobian_policy_t> jacPolicy_;
};

}}}}
#endif  // ODE_IMPLICIT_IMPL_ODE_IMPLICIT_STEPPER_ARBITRARY_IMPL_HPP_
