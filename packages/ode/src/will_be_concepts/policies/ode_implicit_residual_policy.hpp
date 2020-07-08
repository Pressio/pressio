
#ifndef ODE_IMPLICIT_POLICIES_IMPLICIT_RESIDUAL_POLICY_HPP_
#define ODE_IMPLICIT_POLICIES_IMPLICIT_RESIDUAL_POLICY_HPP_

namespace pressio{ namespace ode{ namespace concepts {

template<
  typename T,
  typename tag,
  std::size_t numPrevStates,
  typename state_t,
  typename residual_t,
  typename system_t,
  typename scalar_t,
  typename enable = void
  >
struct implicit_residual_policy : std::false_type{};


template<
  typename T,
  typename tag,
  std::size_t numPrevStates,
  typename state_t,
  typename residual_t,
  typename system_t,
  typename scalar_t
  >
struct implicit_residual_policy<
  T, tag, numPrevStates, state_t, residual_t, system_t, scalar_t,
  ::pressio::mpl::enable_if_t<
    // is callable with two args
    std::is_same<
      residual_t,
      decltype
      (
       std::declval<T const>().create(std::declval<system_t const &>())
       )
      >::value
    and

    // is callable with six
    std::is_void<
      decltype
      (
       std::declval<T const>().template compute
       <tag>(
	     std::declval<state_t const &>(),
	     std::declval<::pressio::ode::AuxStatesManager<state_t, numPrevStates> const &>(),
	     std::declval<system_t const &>(),
	     std::declval<scalar_t const &>(),
	     std::declval<scalar_t const &>(),
	     std::declval<::pressio::ode::types::step_t>(),
	     std::declval<residual_t &>(),
	     ::pressio::Norm::Undefined,
	     std::declval<scalar_t &>()
	     )
       )
      >::value
    >
  > : std::true_type{};
//------------------------------------------------------------------


template<typename T, typename ... args>
using implicit_euler_residual_policy =
  implicit_residual_policy<
  T, ::pressio::ode::implicitmethods::Euler, 1, args...>;

template<typename T, typename ... args>
using implicit_bdf2_residual_policy =
  implicit_residual_policy<
  T, ::pressio::ode::implicitmethods::BDF2, 2, args...>;
//------------------------------------------------------------------

}}} // namespace pressio::ode::concepts
#endif
