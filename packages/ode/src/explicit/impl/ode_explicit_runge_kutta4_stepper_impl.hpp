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

#ifndef ODE_STEPPERS_EXPLICIT_STEPPERS_IMPL_EXPLICIT_RUNGEKUTTA4_STEPPER_IMPL_HPP_
#define ODE_STEPPERS_EXPLICIT_STEPPERS_IMPL_EXPLICIT_RUNGEKUTTA4_STEPPER_IMPL_HPP_

namespace pressio{ namespace ode{ namespace explicitmethods{ namespace impl{

template<
  typename scalar_type,
  typename state_type,
  typename system_type,
  typename velocity_type,
  typename velocity_policy_type,
  typename standard_velocity_policy_type,
  typename ops_t
  >
class ExplicitRungeKutta4StepperImpl
  : public StepperBase< 
  ExplicitRungeKutta4StepperImpl<scalar_type, state_type, system_type, 
    velocity_type, velocity_policy_type, standard_velocity_policy_type, ops_t> 
  >
{

  using this_t = ExplicitRungeKutta4StepperImpl<scalar_type, state_type, system_type, 
    velocity_type, velocity_policy_type, standard_velocity_policy_type, ops_t>;
  using base_t    = StepperBase<this_t>;
  // need to friend base to allow it to access the private methods below
  friend base_t;

public:
  static constexpr bool is_implicit = false;
  static constexpr bool is_explicit = true;

  static constexpr types::stepper_order_t order_value = 4;

private:
  using state_storage_t	    = ::pressio::ode::AuxStatesContainer<is_explicit, state_type, 1>;
  using velocity_storage_t  = VelocitiesContainer<velocity_type, 4>;
  using system_wrapper_t    = ::pressio::ode::impl::OdeSystemWrapper<system_type>;

  state_storage_t stateAuxStorage_;
  velocity_storage_t veloAuxStorage_;
  system_wrapper_t sys_;

  typename std::conditional<
    std::is_same<velocity_policy_type, standard_velocity_policy_type>::value,
    const standard_velocity_policy_type, 
    const velocity_policy_type & 
  >::type policy_;

  const ops_t * udOps_ = nullptr;

public:
  ExplicitRungeKutta4StepperImpl() = delete;
  ~ExplicitRungeKutta4StepperImpl() = default;
  ExplicitRungeKutta4StepperImpl(const ExplicitRungeKutta4StepperImpl & other)  = delete;
  ExplicitRungeKutta4StepperImpl & operator=(const ExplicitRungeKutta4StepperImpl & other)  = delete;
  ExplicitRungeKutta4StepperImpl(ExplicitRungeKutta4StepperImpl && other)  = delete;
  ExplicitRungeKutta4StepperImpl & operator=(ExplicitRungeKutta4StepperImpl && other)  = delete;

  // the following cnstr is enabled if we are using pressio ops
  template <
    typename _ops_t = ops_t, 
    mpl::enable_if_t< std::is_void<_ops_t>::value, int > = 0
  >
  ExplicitRungeKutta4StepperImpl(const state_type & state, 
                                 const system_type & model, 
                                 const velocity_policy_type & policy)
    : stateAuxStorage_{state}, 
      sys_(model), 
      policy_(policy), 
      veloAuxStorage_(policy_.create(model))
  {}  

  // the following cnstr is enabled if we are using user-defined ops
  template <
    typename _ops_t = ops_t, 
    mpl::enable_if_t< !std::is_void<_ops_t>::value, int > = 0
  >
  ExplicitRungeKutta4StepperImpl(const state_type & state, 
                                 const system_type & model, 
                                 const velocity_policy_type & policy,
                                 const _ops_t & udOps)
    : stateAuxStorage_{state},
      sys_(model), 
      policy_(policy), 
      veloAuxStorage_(policy_.create(model)),
      udOps_(&udOps)
  {}  

  // the following cnstr is only enabled if 
  // policy is default constructible and we are using pressio ops
  template <
    typename _vel_pol_t = velocity_policy_type,
    typename _ops_t = ops_t, 
    mpl::enable_if_t< 
      std::is_void<_ops_t>::value and 
      std::is_default_constructible<_vel_pol_t>::value, 
      int > = 0
  >
  ExplicitRungeKutta4StepperImpl(const state_type & state, 
                                 const system_type & model)
    : stateAuxStorage_{state},
      sys_(model), 
      policy_(), // default construct policy
      veloAuxStorage_(policy_.create(model))
  {}

  // the following cnstr is only enabled if 
  // policy is default constructible and we are using user-defined ops
  template <
    typename _vel_pol_t = velocity_policy_type,
    typename _ops_t = ops_t, 
    mpl::enable_if_t< 
      !std::is_void<_ops_t>::value and 
      std::is_default_constructible<_vel_pol_t>::value, 
      int > = 0
  >
  ExplicitRungeKutta4StepperImpl(const state_type & state, 
                                 const system_type & model, 
                                 const _ops_t & udOps)
    : stateAuxStorage_{state},
      sys_(model), 
      policy_(), // default construct policy
      veloAuxStorage_(policy_.create(model)),
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
    auto & ytmp	   = stateAuxStorage_(0);
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
    policy_.compute(stateInOut, auxRhs0, sys_.get(), t);
    this->stage_update_impl(ytmp, stateInOut, auxRhs0, dt_half);

    // stage 2: ytmp = y + auxRhs1*dt_half;
    policy_.compute(ytmp, auxRhs1, sys_.get(), t_phalf);
    this->stage_update_impl(ytmp, stateInOut, auxRhs1, dt_half);

    // stage 3: ytmp = y + auxRhs2*dt;
    policy_.compute(ytmp, auxRhs2, sys_.get(), t_phalf);
    this->stage_update_impl(ytmp, stateInOut, auxRhs2, dt);

    // stage 4: y_n += dt/6 * ( k1 + 2 * k2 + 2 * k3 + k4 )
    policy_.compute(ytmp, auxRhs3, sys_.get(), t + dt);
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
    ::pressio::ops::do_update(yIn, stateIn, one, rhsIn, dtValue);
  }

  template<typename rhs_t, typename _ops_t = ops_t>
  mpl::enable_if_t< std::is_void<_ops_t>::value >
  stage_update_impl(state_type & stateIn,
			 const rhs_t & rhsIn0, const rhs_t & rhsIn1,
			 const rhs_t & rhsIn2, const rhs_t & rhsIn3,
			 scalar_type dt6, scalar_type dt3)
  {
    constexpr auto one  = ::pressio::utils::constants<scalar_type>::one();
    ::pressio::ops::do_update(stateIn, one, rhsIn0, dt6, rhsIn1,
			      dt3, rhsIn2, dt3, rhsIn3, dt6);
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
    udOps_->do_update(*yIn.data(), *stateIn.data(), one, *rhsIn.data(), dtValue);
  }

  template<typename rhs_t, typename _ops_t = ops_t>
  mpl::enable_if_t< !std::is_void<_ops_t>::value >
  stage_update_impl(state_type & stateIn,
			 const rhs_t & rhsIn0, const rhs_t & rhsIn1,
			 const rhs_t & rhsIn2, const rhs_t & rhsIn3,
			 scalar_type dt6, scalar_type dt3)
  {
    constexpr auto one  = ::pressio::utils::constants<scalar_type>::one();
    udOps_->do_update(*stateIn.data(), one, *rhsIn0.data(),
		     dt6, *rhsIn1.data(), dt3, *rhsIn2.data(),
		     dt3, *rhsIn3.data(), dt6);
  }
}; //end class

}}}}//end namespace pressio::ode::explicitmethods::impl
#endif
