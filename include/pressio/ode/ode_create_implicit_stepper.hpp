/*
//@HEADER
// ************************************************************************
//
// ode_create_implicit_stepper.hpp
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

#ifndef ODE_ODE_CREATE_IMPLICIT_STEPPER_HPP_
#define ODE_ODE_CREATE_IMPLICIT_STEPPER_HPP_


#include "./impl/ode_implicit_discrete_residual.hpp"
#include "./impl/ode_implicit_discrete_jacobian.hpp"
#include "./impl/ode_implicit_policy_residual_jacobian_without_mass_matrix.hpp"
#include "./impl/ode_implicit_policy_residual_jacobian_with_mass_matrix.hpp"
#include "./impl/ode_implicit_stepper_standard.hpp"
#include "./impl/ode_implicit_stepper_arbitrary.hpp"
#include "./impl/ode_implicit_create_impl.hpp"

namespace pressio{ namespace ode{

#if defined PRESSIO_ENABLE_CXX20
template<class SystemType>
  requires RealValuedOdeSystemFusingRhsAndJacobian<mpl::remove_cvref_t<SystemType>>
  && (Traits<typename mpl::remove_cvref_t<SystemType>::state_type>::rank    == 1)
  && (Traits<typename mpl::remove_cvref_t<SystemType>::rhs_type>::rank      == 1)
  && (Traits<typename mpl::remove_cvref_t<SystemType>::jacobian_type>::rank == 2)
  && requires(      typename mpl::remove_cvref_t<SystemType>::state_type    & s,
	            typename mpl::remove_cvref_t<SystemType>::rhs_type      & r,
	            typename mpl::remove_cvref_t<SystemType>::jacobian_type & J,
	      const typename mpl::remove_cvref_t<SystemType>::state_type & s1,
	      const typename mpl::remove_cvref_t<SystemType>::state_type & s2,
	      const typename mpl::remove_cvref_t<SystemType>::state_type & s3,
	      const typename mpl::remove_cvref_t<SystemType>::rhs_type   & r1,
	      const typename mpl::remove_cvref_t<SystemType>::rhs_type   & r2,
	      ode::scalar_of_t< mpl::remove_cvref_t<SystemType> > a,
	      ode::scalar_of_t< mpl::remove_cvref_t<SystemType> > b,
	      ode::scalar_of_t< mpl::remove_cvref_t<SystemType> > c,
  	      ode::scalar_of_t< mpl::remove_cvref_t<SystemType> > d,
	      ode::scalar_of_t< mpl::remove_cvref_t<SystemType> > e)
  {
    { ::pressio::ops::deep_copy(s, s1) };

    // bdf1, bdf2
    { ::pressio::ops::update(r, a, s1, b, s2, c) };
    { ::pressio::ops::update(r, a, s1, b, s2, c, s3, d) };
    // cn
    { ::pressio::ops::update(r, a, s1, b, s2, c, r1, d, r2, e) };

    // all
    { ::pressio::ops::scale(J, a) };
    { ::pressio::ops::add_to_diagonal(J, a) };
  }
#else
template<
  class SystemType,
  mpl::enable_if_t<
    RealValuedOdeSystemFusingRhsAndJacobian<mpl::remove_cvref_t<SystemType>>::value,
    int > = 0
  >
#endif
auto create_implicit_stepper(StepScheme schemeName,                     // (1)
			     SystemType && system)
{

  assert(schemeName == StepScheme::BDF1 ||
	 schemeName == StepScheme::BDF2 ||
	 schemeName == StepScheme::CrankNicolson);

  using system_type   = mpl::remove_cvref_t<SystemType>;
  using ind_var_type  = typename system_type::independent_variable_type;
  using state_type    = typename system_type::state_type;
  using residual_type = typename system_type::rhs_type;
  using jacobian_type = typename system_type::jacobian_type;

  // it is very important to use "SystemType" as template arg
  // because that it the right type carrying how we store the system
  using policy_type = impl::ResidualJacobianStandardPolicy<
    SystemType, ind_var_type, state_type, residual_type, jacobian_type>;

  using impl_type = impl::ImplicitStepperStandardImpl<
    ind_var_type, state_type, residual_type,
    jacobian_type, policy_type>;
  return impl::create_implicit_stepper_impl<
    impl_type>(schemeName, policy_type(std::forward<SystemType>(system)));
}



#if defined PRESSIO_ENABLE_CXX20
template<class SystemType>
  requires RealValuedCompleteOdeSystem<mpl::remove_cvref_t<SystemType>>
  && (Traits<typename mpl::remove_cvref_t<SystemType>::state_type>::rank       == 1)
  && (Traits<typename mpl::remove_cvref_t<SystemType>::rhs_type>::rank         == 1)
  && (Traits<typename mpl::remove_cvref_t<SystemType>::jacobian_type>::rank    == 2)
  && (Traits<typename mpl::remove_cvref_t<SystemType>::mass_matrix_type>::rank == 2)
  && requires(      typename mpl::remove_cvref_t<SystemType>::state_type       & s,
	            typename mpl::remove_cvref_t<SystemType>::rhs_type         & r,
	            typename mpl::remove_cvref_t<SystemType>::jacobian_type    & J,
	      const typename mpl::remove_cvref_t<SystemType>::mass_matrix_type & M,
	      const typename mpl::remove_cvref_t<SystemType>::state_type & s1,
	      const typename mpl::remove_cvref_t<SystemType>::state_type & s2,
	      const typename mpl::remove_cvref_t<SystemType>::state_type & s3,
	      const typename mpl::remove_cvref_t<SystemType>::rhs_type & r1,
	      ode::scalar_of_t< mpl::remove_cvref_t<SystemType> > a,
	      ode::scalar_of_t< mpl::remove_cvref_t<SystemType> > b,
	      ode::scalar_of_t< mpl::remove_cvref_t<SystemType> > c,
	      ode::scalar_of_t< mpl::remove_cvref_t<SystemType> > d)
  {
    { ::pressio::ops::deep_copy(s, s1) };

    // bdf1 and bdf2
    { ::pressio::ops::update(s, a, s1, b, s2, c) };
    { ::pressio::ops::update(s, a, s1, b, s2, c, s3, d) };
    { ::pressio::ops::product(::pressio::nontranspose(), a, M, s1, b, r) };
    { ::pressio::ops::update(r, a, r1, b) };
    { ::pressio::ops::update(J, a, M, b)  };
  }
#else
template<
  class SystemType,
  mpl::enable_if_t<
    RealValuedCompleteOdeSystem<mpl::remove_cvref_t<SystemType>>::value,
    int > = 0
  >
#endif
auto create_implicit_stepper(StepScheme schemeName,                     // (2)
			     SystemType && system)
{

  assert(schemeName == StepScheme::BDF1 ||
	 schemeName == StepScheme::BDF2);

  using system_type   = mpl::remove_cvref_t<SystemType>;
  using ind_var_type  = typename system_type::independent_variable_type;
  using state_type    = typename system_type::state_type;
  using residual_type = typename system_type::rhs_type;
  using jacobian_type = typename system_type::jacobian_type;
  using mass_mat_type = typename system_type::mass_matrix_type;

  // it is very important to use "SystemType" as template arg
  // because that it the right type carrying how we store the system
  using policy_type = impl::ResidualJacobianWithMassMatrixStandardPolicy<
    SystemType, ind_var_type, state_type,
    residual_type, jacobian_type, mass_mat_type>;

  using impl_type = impl::ImplicitStepperStandardImpl<
    ind_var_type, state_type, residual_type,
    jacobian_type, policy_type>;
  return impl::create_implicit_stepper_impl<
    impl_type>(schemeName, policy_type(std::forward<SystemType>(system)));
}


template<
  class ResidualJacobianPolicyType
#if not defined PRESSIO_ENABLE_CXX20
  ,mpl::enable_if_t<
    ::pressio::ode::ImplicitResidualJacobianPolicy<
      mpl::remove_cvref_t<ResidualJacobianPolicyType>>::value, int
    > = 0
#endif
  >
#ifdef PRESSIO_ENABLE_CXX20
requires ::pressio::ode::ImplicitResidualJacobianPolicy<
  mpl::remove_cvref_t<ResidualJacobianPolicyType>
  >
#endif
auto create_implicit_stepper(StepScheme schemeName,
			     ResidualJacobianPolicyType && policy)
{

  assert(schemeName == StepScheme::BDF1 ||
	 schemeName == StepScheme::BDF2 ||
	 schemeName == StepScheme::CrankNicolson);

  using policy_type   = mpl::remove_cvref_t<ResidualJacobianPolicyType>;
  using ind_var_type  = typename policy_type::independent_variable_type;
  using state_type    = typename policy_type::state_type;
  using residual_type = typename policy_type::residual_type;
  using jacobian_type = typename policy_type::jacobian_type;

  using impl_type = impl::ImplicitStepperStandardImpl<
    ind_var_type, state_type, residual_type,
    jacobian_type, ResidualJacobianPolicyType>;

  return impl::create_implicit_stepper_impl<
    impl_type>(schemeName, std::forward<ResidualJacobianPolicyType>(policy));
}

//
// auxiliary API
//
template<class ...Args>
auto create_bdf1_stepper(Args && ... args){
  return create_implicit_stepper(StepScheme::BDF1,
				 std::forward<Args>(args)...);
}

template<class ...Args>
auto create_bdf2_stepper(Args && ... args){
  return create_implicit_stepper(StepScheme::BDF2,
				 std::forward<Args>(args)...);
}

template<class ...Args>
auto create_cranknicolson_stepper(Args && ... args){
  return create_implicit_stepper(StepScheme::CrankNicolson,
				 std::forward<Args>(args)...);
}


//
// num of states as template arg constructs the arbitrary stepper
//
template<int TotalNumberOfDesiredStates, class SystemType>
#if defined PRESSIO_ENABLE_CXX20
  requires RealValuedFullyDiscreteSystemWithJacobian<
     mpl::remove_cvref_t<SystemType>, TotalNumberOfDesiredStates>
  && (Traits<typename mpl::remove_cvref_t<SystemType>::state_type>::rank == 1)
  && (Traits<typename mpl::remove_cvref_t<SystemType>::discrete_residual_type>::rank == 1)
  && (Traits<typename mpl::remove_cvref_t<SystemType>::discrete_jacobian_type>::rank == 2)
  && requires(      typename mpl::remove_cvref_t<SystemType>::state_type & s,
	            typename mpl::remove_cvref_t<SystemType>::discrete_residual_type & r,
	            typename mpl::remove_cvref_t<SystemType>::discrete_jacobian_type & J,
	      const typename mpl::remove_cvref_t<SystemType>::state_type & s1,
	      const typename mpl::remove_cvref_t<SystemType>::state_type & s2,
	      const typename mpl::remove_cvref_t<SystemType>::state_type & s3,
	      ode::scalar_of_t< mpl::remove_cvref_t<SystemType>, TotalNumberOfDesiredStates > a)
  {
    { ::pressio::ops::deep_copy(s, s1) };
  }
#endif
auto create_implicit_stepper(SystemType && system)
{

  using system_type = mpl::remove_cvref_t<SystemType>;
#if not defined PRESSIO_ENABLE_CXX20
  static_assert(RealValuedFullyDiscreteSystemWithJacobian<system_type, TotalNumberOfDesiredStates>::value,
		"The system passed does not meet the FullyDiscrete API");
#endif

  using ind_var_type  = typename system_type::independent_variable_type;
  using state_type    = typename system_type::state_type;
  using residual_type = typename system_type::discrete_residual_type;
  using jacobian_type = typename system_type::discrete_jacobian_type;

  using stepper_type = impl::StepperArbitrary<
    TotalNumberOfDesiredStates, SystemType, ind_var_type,
    state_type, residual_type, jacobian_type
    >;
  return stepper_type(std::forward<SystemType>(system));
}

}} // end namespace pressio::ode
#endif  // ODE_ODE_CREATE_IMPLICIT_STEPPER_HPP_
