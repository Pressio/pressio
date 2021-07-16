/*
//@HEADER
// ************************************************************************
//
// ode_explicit_adams_bashforth2_stepper_impl.hpp
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

#ifndef ODE_EXPLICIT_IMPL_ODE_EXPLICIT_ADAMS_BASHFORTH2_STEPPER_IMPL_HPP_
#define ODE_EXPLICIT_IMPL_ODE_EXPLICIT_ADAMS_BASHFORTH2_STEPPER_IMPL_HPP_

namespace pressio{ namespace ode{ namespace explicitmethods{ namespace impl{

template<
  typename scalar_type,
  typename state_type,
  typename system_type,
  typename velocity_type,
  typename velocity_policy_type,
  typename ops_t,
  bool is_standard_policy
  >
class ExplicitAdamsBashforth2Stepper
{
public:
  static constexpr bool is_implicit = false;
  static constexpr bool is_explicit = true;
  static constexpr types::stepper_order_t order_value = 2;
  using velocity_storage_t  = ::pressio::containers::IndexableStaticCollection<velocity_type, 2>;

private:
  std::reference_wrapper<const system_type> systemObj_;
  ::pressio::utils::instance_or_reference_wrapper<velocity_policy_type> policy_;
  velocity_storage_t velocities_;
  state_type tmpState_;
  const ops_t * udOps_ = nullptr;

public:
  ExplicitAdamsBashforth2Stepper() = delete;
  ExplicitAdamsBashforth2Stepper(const ExplicitAdamsBashforth2Stepper & other) = default;
  ExplicitAdamsBashforth2Stepper & operator=(const ExplicitAdamsBashforth2Stepper & other) = delete;
  ExplicitAdamsBashforth2Stepper(ExplicitAdamsBashforth2Stepper && other)  = default;
  ExplicitAdamsBashforth2Stepper & operator=(ExplicitAdamsBashforth2Stepper && other)  = delete;
  ~ExplicitAdamsBashforth2Stepper() = default;

  template <
    typename _ops_t = ops_t,
    mpl::enable_if_t< std::is_void<_ops_t>::value, int > = 0
    >
  ExplicitAdamsBashforth2Stepper(const state_type & state,
				 const system_type & systemObj,
				 const mpl::remove_cvref_t<velocity_policy_type> & policy)
    : systemObj_(systemObj),
      policy_(policy),
      velocities_(policy_.get().create(systemObj)),
      tmpState_{state}
  {}

  template <
    typename _ops_t = ops_t,
    mpl::enable_if_t< !std::is_void<_ops_t>::value, int > = 0
    >
  ExplicitAdamsBashforth2Stepper(const state_type & state,
				 const system_type & systemObj,
				 const mpl::remove_cvref_t<velocity_policy_type> & policy,
				 const _ops_t & udOps)
    : systemObj_(systemObj),
      policy_(policy),
      velocities_(policy_.get().create(systemObj)),
      tmpState_{state},
      udOps_(&udOps)
  {}

  // only enabled if policy standard and using pressio ops
  template <
    bool _is_standard_policy = is_standard_policy,
    typename _ops_t = ops_t,
    mpl::enable_if_t<
      _is_standard_policy and std::is_void<_ops_t>::value,
      int > = 0
    >
  ExplicitAdamsBashforth2Stepper(const state_type & state,
				 const system_type & systemObj)
    : systemObj_(systemObj),
      policy_(),
      velocities_(policy_.get().create(systemObj)),
      tmpState_{state}
  {}

  // only enabled if policy standard and user-defined ops
  template <
    bool _is_standard_policy = is_standard_policy,
    typename _ops_t = ops_t,
    mpl::enable_if_t<
      _is_standard_policy and !std::is_void<_ops_t>::value,
      int > = 0
    >
  ExplicitAdamsBashforth2Stepper(const state_type & state,
				 const system_type & systemObj,
				 const _ops_t & udOps)
    : systemObj_(systemObj),
      policy_(),
      velocities_(policy_.get().create(systemObj)),
      tmpState_{state},
      udOps_(&udOps)
  {}

public:
  types::stepper_order_t order() const
  {
    return order_value;
  }

  void doStep(state_type & odeSolution,
  	      const scalar_type & currentTime,
  	      const scalar_type & dt,
  	      const types::step_t & stepNumber)
  {
    PRESSIOLOG_DEBUG("adams-bashforth2 stepper: do step");

    // y_n+1 = y_n + dt*[ (3/2)*f(y_n, t_n) - (1/2)*f(y_n-1, t_n-1) ]

    const auto cfn   = ::pressio::utils::constants<scalar_type>::threeOvTwo()*dt;
    const auto cfnm1 = ::pressio::utils::constants<scalar_type>::negOneHalf()*dt;

    if (stepNumber==1){
      // use Euler forward or we could use something else here maybe RK4
      auto & rhs = velocities_(0);
      policy_.get().compute(odeSolution, rhs, systemObj_.get(), currentTime);
      doUpdate1(odeSolution, rhs, dt);
    }
    else{
      auto & fn   = velocities_(0);
      auto & fnm1 = velocities_(1);
      // fn -> fnm1
      updateStoredVelocities(fnm1, fn);

      policy_.get().compute(odeSolution, fn, systemObj_.get(), currentTime);
      doUpdate2(odeSolution, fn, cfn, fnm1, cfnm1);
    }
  }

private:
  template<typename T, typename _ops_t = ops_t>
  mpl::enable_if_t< std::is_void<_ops_t>::value >
  updateStoredVelocities(T & to, const T & from)
  {
    ::pressio::ops::deep_copy(to, from);
  }

  template<typename T, typename _ops_t = ops_t>
  mpl::enable_if_t< mpl::not_void<_ops_t>::value >
  updateStoredVelocities(T & to, const T & from)
  {
    udOps_->deep_copy(*to.data(), *from.data());
  }

  template<typename f_t, typename _ops_t = ops_t>
  mpl::enable_if_t< std::is_void<_ops_t>::value >
  doUpdate1(state_type & odeSolution,
	    const f_t & rhs,
	    const scalar_type & dt)
  {
    constexpr auto one   = ::pressio::utils::constants<scalar_type>::one();
    ::pressio::ops::update(odeSolution, one, rhs, dt);
  }

  template<typename f_t, typename _ops_t = ops_t>
  mpl::enable_if_t< mpl::not_void<_ops_t>::value >
  doUpdate1(state_type & odeSolution,
	    const f_t & rhs,
	    const scalar_type & dt)
  {
    constexpr auto one   = ::pressio::utils::constants<scalar_type>::one();
    udOps_->update(*odeSolution.data(), one, *rhs.data(), dt);
  }

  template<typename f_t, typename _ops_t = ops_t>
  mpl::enable_if_t< std::is_void<_ops_t>::value >
  doUpdate2(state_type & odeSolution,
	    const f_t & fn, const scalar_type & cfn,
	    const f_t & fnm1, const scalar_type & cfnm1)
  {
    constexpr auto one   = ::pressio::utils::constants<scalar_type>::one();
    ::pressio::ops::update(odeSolution, one, fn, cfn, fnm1, cfnm1);
  }

  template<typename f_t, typename _ops_t = ops_t>
  mpl::enable_if_t< mpl::not_void<_ops_t>::value >
  doUpdate2(state_type & odeSolution,
	    const f_t & fn, const scalar_type & cfn,
	    const f_t & fnm1, const scalar_type & cfnm1)
  {
    constexpr auto one   = ::pressio::utils::constants<scalar_type>::one();
    udOps_->update(*odeSolution.data(), one,
		   *fn.data(), cfn,
		   *fnm1.data(), cfnm1);
  }
};

}}}}//end namespace pressio::ode::explicitmethods::impl
#endif  // ODE_EXPLICIT_IMPL_ODE_EXPLICIT_ADAMS_BASHFORTH2_STEPPER_IMPL_HPP_
