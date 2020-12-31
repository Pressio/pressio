/*
//@HEADER
// ************************************************************************
//
// ode_explicit_runge_kutta4_stepper_impl.hpp
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

#ifndef ODE_EXPLICIT_IMPL_ODE_EXPLICIT_RUNGE_KUTTA4_STEPPER_IMPL_HPP_
#define ODE_EXPLICIT_IMPL_ODE_EXPLICIT_RUNGE_KUTTA4_STEPPER_IMPL_HPP_

namespace pressio{ namespace ode{ namespace explicitmethods{ namespace impl{

template<
  typename scalar_type,
  typename state_type,
  typename system_type,
  typename velocity_type,
  typename velocity_policy_type,
  typename ops_t
  >
class ExplicitRungeKutta4Stepper
{
public:
  static constexpr bool is_implicit = false;
  static constexpr bool is_explicit = true;
  static constexpr types::stepper_order_t order_value = 4;
  using velocity_storage_t  = ::pressio::containers::IndexableStaticCollection<velocity_type, 4>;

  static constexpr bool using_standard_policy =
    std::is_same<
    velocity_policy_type,
    ::pressio::ode::explicitmethods::policy::VelocityStandardPolicy<state_type>
    >::value;

private:
  std::reference_wrapper<const system_type> systemObj_;
  ::pressio::utils::instance_or_reference_wrapper<velocity_policy_type> policy_;
  velocity_storage_t veloAuxStorage_;
  state_type tmpState_;
  const ops_t * udOps_ = nullptr;

public:
  ExplicitRungeKutta4Stepper() = delete;
  ExplicitRungeKutta4Stepper(const ExplicitRungeKutta4Stepper & other) = default;
  ExplicitRungeKutta4Stepper & operator=(const ExplicitRungeKutta4Stepper & other) = delete;
  ExplicitRungeKutta4Stepper(ExplicitRungeKutta4Stepper && other)  = default;
  ExplicitRungeKutta4Stepper & operator=(ExplicitRungeKutta4Stepper && other)  = delete;
  ~ExplicitRungeKutta4Stepper() = default;

  // the following cnstr is enabled if we are using pressio ops
  // cnstr enabled if we are using pressio ops
  template <
    typename _ops_t = ops_t,
    mpl::enable_if_t< std::is_void<_ops_t>::value, int > = 0
    >
  ExplicitRungeKutta4Stepper(const state_type & state,
			     const system_type & systemObj,
			     const mpl::remove_cvref_t<velocity_policy_type> & policy)
    : systemObj_(systemObj),
      policy_(policy),
      veloAuxStorage_(policy_.get().create(systemObj)),
      tmpState_{state}
  {}

  // the following cnstr is enabled if we are using user-defined ops
  template <
    typename _ops_t = ops_t,
    mpl::enable_if_t< !std::is_void<_ops_t>::value, int > = 0
    >
  ExplicitRungeKutta4Stepper(const state_type & state,
			     const system_type & systemObj,
			     const mpl::remove_cvref_t<velocity_policy_type> & policy,
			     const _ops_t & udOps)
    : systemObj_(systemObj),
      policy_(policy),
      veloAuxStorage_(policy_.get().create(systemObj)),
      tmpState_{state},
      udOps_(&udOps)
  {}

  // only enabled if policy standard and using pressio ops
  template <
    bool _using_standard_policy = using_standard_policy,
    typename _ops_t = ops_t,
    mpl::enable_if_t<
      _using_standard_policy and std::is_void<_ops_t>::value,
      int > = 0
    >
  ExplicitRungeKutta4Stepper(const state_type & state,
			     const system_type & systemObj)
    : systemObj_(systemObj),
      policy_(),
      veloAuxStorage_(policy_.get().create(systemObj)),
      tmpState_{state}
  {}

  // only enabled if policy standard and user-defined ops
  template <
    bool _using_standard_policy = using_standard_policy,
    typename _ops_t = ops_t,
    mpl::enable_if_t<
      _using_standard_policy and !std::is_void<_ops_t>::value,
      int > = 0
    >
  ExplicitRungeKutta4Stepper(const state_type & state,
			     const system_type & systemObj,
			     const _ops_t & udOps)
    : systemObj_(systemObj),
      policy_(),
      veloAuxStorage_(policy_.get().create(systemObj)),
      tmpState_{state},
      udOps_(&udOps)
  {}

public:
  types::stepper_order_t order() const
  {
    return order_value;
  }

  void doStep(state_type & stateInOut,
  	      const scalar_type & t,
  	      const scalar_type & dt,
  	      const types::step_t & step)
  {
    PRESSIOLOG_DEBUG("rk4 stepper: do step");

    auto & auxRhs0 = veloAuxStorage_(0);
    auto & auxRhs1 = veloAuxStorage_(1);
    auto & auxRhs2 = veloAuxStorage_(2);
    auto & auxRhs3 = veloAuxStorage_(3);

    constexpr auto two  = ::pressio::utils::constants<scalar_type>::two();
    constexpr auto three  = ::pressio::utils::constants<scalar_type>::three();
    constexpr auto six  = two * three;

    const scalar_type dt_half = dt / two;
    const scalar_type t_phalf = t + dt_half;
    const scalar_type dt6 = dt / six;
    const scalar_type dt3 = dt / three;

    // stage 1: ytmp = y + auxRhs0*dt_half;
    policy_.get().compute(stateInOut, auxRhs0, systemObj_.get(), t);
    this->stage_update_impl(tmpState_, stateInOut, auxRhs0, dt_half);

    // stage 2: ytmp = y + auxRhs1*dt_half;
    policy_.get().compute(tmpState_, auxRhs1, systemObj_.get(), t_phalf);
    this->stage_update_impl(tmpState_, stateInOut, auxRhs1, dt_half);

    // stage 3: ytmp = y + auxRhs2*dt;
    policy_.get().compute(tmpState_, auxRhs2, systemObj_.get(), t_phalf);
    this->stage_update_impl(tmpState_, stateInOut, auxRhs2, dt);

    // stage 4: y_n += dt/6 * ( k1 + 2 * k2 + 2 * k3 + k4 )
    policy_.get().compute(tmpState_, auxRhs3, systemObj_.get(), t + dt);
    this->stage_update_impl(stateInOut, auxRhs0, auxRhs1, auxRhs2, auxRhs3, dt6, dt3);
  }//end doStep

private:
  /* use pressio ops  */
  template<typename rhs_t, typename _ops_t = ops_t>
  mpl::enable_if_t<std::is_void<_ops_t>::value>
  stage_update_impl(state_type & yIn,
		    const state_type & stateIn,
		    const rhs_t & rhsIn,
		    scalar_type dtValue)
  {
    constexpr auto one  = ::pressio::utils::constants<scalar_type>::one();
    ::pressio::ops::update(yIn, stateIn, one, rhsIn, dtValue);
  }

  template<typename rhs_t, typename _ops_t = ops_t>
  mpl::enable_if_t< std::is_void<_ops_t>::value >
  stage_update_impl(state_type & stateIn,
		    const rhs_t & rhsIn0, const rhs_t & rhsIn1,
		    const rhs_t & rhsIn2, const rhs_t & rhsIn3,
		    scalar_type dt6, scalar_type dt3)
  {
    constexpr auto one  = ::pressio::utils::constants<scalar_type>::one();
    ::pressio::ops::update(stateIn, one,
			   rhsIn0, dt6,
			   rhsIn1, dt3,
			   rhsIn2, dt3,
			   rhsIn3, dt6);
  }
  // -------------------------------------------------------

  /* with user defined ops */
  template<typename rhs_t, typename _ops_t = ops_t>
  mpl::enable_if_t< !std::is_void<_ops_t>::value >
  stage_update_impl(state_type & yIn,
		    const state_type & stateIn,
		    const rhs_t & rhsIn,
		    scalar_type dtValue)
  {
    constexpr auto one  = ::pressio::utils::constants<scalar_type>::one();
    udOps_->update(*yIn.data(), *stateIn.data(), one, *rhsIn.data(), dtValue);
  }

  template<typename rhs_t, typename _ops_t = ops_t>
  mpl::enable_if_t< !std::is_void<_ops_t>::value >
  stage_update_impl(state_type & stateIn,
		    const rhs_t & rhsIn0, const rhs_t & rhsIn1,
		    const rhs_t & rhsIn2, const rhs_t & rhsIn3,
		    scalar_type dt6, scalar_type dt3)
  {
    constexpr auto one  = ::pressio::utils::constants<scalar_type>::one();
    udOps_->update(*stateIn.data(), one, *rhsIn0.data(),
		   dt6, *rhsIn1.data(), dt3, *rhsIn2.data(),
		   dt3, *rhsIn3.data(), dt6);
  }
};

}}}}//end namespace pressio::ode::explicitmethods::impl
#endif  // ODE_EXPLICIT_IMPL_ODE_EXPLICIT_RUNGE_KUTTA4_STEPPER_IMPL_HPP_
