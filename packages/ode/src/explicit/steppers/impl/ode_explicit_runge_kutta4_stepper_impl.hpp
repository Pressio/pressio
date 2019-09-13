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

#include "../ode_explicit_stepper_base.hpp"

namespace pressio{ namespace ode{ namespace impl{

template<
  typename scalar_type,
  typename state_type,
  typename system_type,
  typename ode_residual_type,
  typename residual_policy_type,
  typename ops_t
  >
class ExplicitRungeKutta4StepperImpl<scalar_type,
				     state_type,
				     system_type,
				     ode_residual_type,
				     residual_policy_type,
				     ops_t>
{

  static_assert( meta::is_legitimate_explicit_velocity_policy<
		 residual_policy_type>::value ||
		 meta::is_explicit_runge_kutta4_velocity_standard_policy<
		 residual_policy_type>::value,
"EXPLICIT RUNGEKUTTA4 VELOCITY_POLICY NOT ADMISSIBLE, \
MAYBE NOT A CHILD OF ITS BASE OR DERIVING FROM WRONG BASE");

  using this_t = ExplicitRungeKutta4StepperImpl< scalar_type,
						 state_type, system_type,
						 ode_residual_type,
						 residual_policy_type,
						 ops_t>;

  using state_storage_t = OdeStorage<state_type, 1>;
  using resid_storage_t = OdeStorage<ode_residual_type, 4>;
  using system_wrapper_t = OdeSystemWrapper<system_type>;

  state_storage_t stateAuxStorage_;
  resid_storage_t residAuxStorage_;
  system_wrapper_t sys_;
  const residual_policy_type & policy_;

public:
  ExplicitRungeKutta4StepperImpl(const system_type & model,
  				 const residual_policy_type & res_policy_obj,
  				 const state_type & y0,
  				 const ode_residual_type & r0)
    : stateAuxStorage_{y0}, residAuxStorage_{r0},
      sys_{model}, policy_{res_policy_obj}{}

  ExplicitRungeKutta4StepperImpl() = delete;
  ~ExplicitRungeKutta4StepperImpl() = default;

public:
  template< typename step_t >
  void doStep(state_type & y,
  	      scalar_type t,
  	      scalar_type dt,
  	      step_t step)
  {
    auto & ytmp	   = stateAuxStorage_.data_[0];
    auto & auxRhs0 = residAuxStorage_.data_[0];
    auto & auxRhs1 = residAuxStorage_.data_[1];
    auto & auxRhs2 = residAuxStorage_.data_[2];
    auto & auxRhs3 = residAuxStorage_.data_[3];

    constexpr auto two  = ::pressio::utils::constants::two<scalar_type>();
    constexpr auto three  = ::pressio::utils::constants::three<scalar_type>();
    constexpr auto six  = two * three;

    const scalar_type dt_half = dt / two;
    const scalar_type t_phalf = t + dt_half;
    const scalar_type dt6 = dt / six;
    const scalar_type dt3 = dt / three;

    // stage 1: ytmp = y + auxRhs0*dt_half;
    policy_(y, auxRhs0, sys_.get(), t);
    this->stage_update_impl(ytmp, y, auxRhs0, dt_half);
    //::pressio::containers::ops::do_update(ytmp, y, one, auxRhs0, dt_half);

    // stage 2: ytmp = y + auxRhs1*dt_half;
    policy_(ytmp, auxRhs1, sys_.get(), t_phalf);
    this->stage_update_impl(ytmp, y, auxRhs1, dt_half);
    //::pressio::containers::ops::do_update(ytmp, y, one, auxRhs1, dt_half);

    // stage 3: ytmp = y + auxRhs2*dt;
    policy_(ytmp, auxRhs2, sys_.get(), t_phalf);
    this->stage_update_impl(ytmp, y, auxRhs2, dt);
    //::pressio::containers::ops::do_update(ytmp, y, one, auxRhs2, dt);

    // stage 4: y_n += dt/6 * ( k1 + 2 * k2 + 2 * k3 + k4 )
    policy_(ytmp, auxRhs3, sys_.get(), t + dt);
    this->stage_update_impl(y, auxRhs0, auxRhs1, auxRhs2, auxRhs3, dt6, dt3);
    //::pressio::containers::ops::do_update(y, one,
					// auxRhs0, dt6,
					// auxRhs1, dt3,
					// auxRhs2, dt3,
					// auxRhs3, dt6);
  }//end doStep

private:
  /* user defined ops are void, use those from containers */
  template<
    typename rhs_t,
    typename _state_type = state_type,
    typename _ops_t = ops_t,
    mpl::enable_if_t<
      std::is_void<_ops_t>::value
      > * = nullptr
  >
  void stage_update_impl(_state_type & yIn,
			 const _state_type & y,
			 const rhs_t & rhsIn,
			 scalar_type dtValue){
    constexpr auto one  = ::pressio::utils::constants::one<scalar_type>();
    ::pressio::containers::ops::do_update(yIn, y, one, rhsIn, dtValue);
  }

  /* with user defined ops */
  template<
    typename rhs_t,
    typename _state_type = state_type,
    typename _ops_t = ops_t,
    mpl::enable_if_t<
      ::pressio::containers::meta::is_wrapper<_state_type>::value and
      !std::is_void<_ops_t>::value
      > * = nullptr
  >
  void stage_update_impl(_state_type & yIn,
			 const _state_type & y,
			 const rhs_t & rhsIn,
			 scalar_type dtValue){
    constexpr auto one  = ::pressio::utils::constants::one<scalar_type>();
    using op = typename ops_t::update_op;
    op::do_update(*yIn.data(), *y.data(), one, *rhsIn.data(), dtValue);
  }

  /* user defined ops are void, use those from containers */
  template<
    typename rhs_t,
    typename _state_type = state_type,
    typename _ops_t = ops_t,
    mpl::enable_if_t<
      std::is_void<_ops_t>::value
      > * = nullptr
  >
  void stage_update_impl(_state_type & y,
			 const rhs_t & rhsIn0,
			 const rhs_t & rhsIn1,
			 const rhs_t & rhsIn2,
			 const rhs_t & rhsIn3,
			 scalar_type dt6,
			 scalar_type dt3){
    constexpr auto one  = ::pressio::utils::constants::one<scalar_type>();
    ::pressio::containers::ops::do_update(y, one, rhsIn0, dt6,
					  rhsIn1, dt3, rhsIn2, dt3,
					  rhsIn3, dt6);
  }

  /* user defined ops */
  template<
    typename rhs_t,
    typename _state_type = state_type,
    typename _ops_t = ops_t,
    mpl::enable_if_t<
      ::pressio::containers::meta::is_wrapper<_state_type>::value and
      !std::is_void<_ops_t>::value
      > * = nullptr
  >
  void stage_update_impl(_state_type & y,
			 const rhs_t & rhsIn0,
			 const rhs_t & rhsIn1,
			 const rhs_t & rhsIn2,
			 const rhs_t & rhsIn3,
			 scalar_type dt6,
			 scalar_type dt3){
    constexpr auto one  = ::pressio::utils::constants::one<scalar_type>();
    using op = typename ops_t::update_op;
    op::do_update(*y.data(), one, *rhsIn0.data(),
		  dt6, *rhsIn1.data(),
		  dt3, *rhsIn2.data(),
		  dt3, *rhsIn3.data(), dt6);
  }
}; //end class

}}}//end namespace pressio::ode::impl
#endif
