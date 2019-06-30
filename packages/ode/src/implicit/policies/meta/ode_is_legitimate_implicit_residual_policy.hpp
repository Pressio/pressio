
#ifndef ODE_IMPLICIT_POLICIES_IS_LEGITIMATE_IMPLICIT_RESIDUAL_POLICY_HPP_
#define ODE_IMPLICIT_POLICIES_IS_LEGITIMATE_IMPLICIT_RESIDUAL_POLICY_HPP_

#include "../base/ode_implicit_residual_policy_base.hpp"

namespace rompp{ namespace ode{ namespace meta {

template<
  typename T,
  ImplicitEnum name,
  int nstates,
  typename state_t,
  typename residual_t,
  typename model_t,
  typename scalar_t,
  typename enable = void
  >
struct is_legitimate_implicit_residual_policy : std::false_type{};


template<
  typename T,
  ImplicitEnum name,
  int nstates,
  typename state_t,
  typename residual_t,
  typename model_t,
  typename scalar_t
  >
struct is_legitimate_implicit_residual_policy<
  T, name, nstates, state_t, residual_t, model_t, scalar_t,
  ::rompp::mpl::enable_if_t<
    // is callable with five args
    std::is_same<
      residual_t,
      decltype
      (
       std::declval<T>().template operator()
       <
       name,
       nstates
       >( std::declval<const state_t &>(),
	  std::declval<const std::array<state_t, nstates> &>(),
	  std::declval<const model_t&>(),
	  std::declval<scalar_t>(),
	  std::declval<scalar_t>()
	  )
       )
      >::value
    and
    // is callable with six
    std::is_void<
      decltype
      (
       std::declval<T>().template operator()
       <
       name,
       nstates
       >( std::declval<const state_t &>(),
	  std::declval<residual_t &>(),
	  std::declval<const std::array<state_t, nstates> &>(),
	  std::declval<const model_t&>(),
	  std::declval<scalar_t>(),
	  std::declval<scalar_t>()
	  )
       )
      >::value
    >
  > : std::true_type{};
//------------------------------------------------------------------


template<typename T, typename ... args>
using is_legitimate_implicit_euler_residual_policy =
  is_legitimate_implicit_residual_policy<
  T, ImplicitEnum::Euler, 1, args...>;

template<typename T, typename ... args>
using is_legitimate_implicit_bdf2_residual_policy =
  is_legitimate_implicit_residual_policy<
  T, ImplicitEnum::BDF2, 2, args...>;
//------------------------------------------------------------------

}}} // namespace rompp::ode::meta
#endif
