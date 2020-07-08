
#ifndef ODE_POLICIES_META_IS_LEGITIMATE_IMPLICIT_JACOBIAN_POLICIES_HPP_
#define ODE_POLICIES_META_IS_LEGITIMATE_IMPLICIT_JACOBIAN_POLICIES_HPP_

namespace pressio{ namespace ode{ namespace concepts {

template<
  typename T,
  typename tag,
  std::size_t numPrevStates,
  typename state_t,
  typename jacobian_t,
  typename system_t,
  typename scalar_t,
  typename = void
  >
struct implicit_jacobian_policy : std::false_type{};

template<
  typename T,
  typename tag,
  std::size_t numPrevStates,
  typename state_t,
  typename jacobian_t,
  typename system_t,
  typename scalar_t
  >
struct implicit_jacobian_policy
<T, tag, numPrevStates, state_t, jacobian_t, system_t, scalar_t,
 ::pressio::mpl::enable_if_t<
   std::is_same<
     jacobian_t,
     decltype
     (
      std::declval<T const>().create(std::declval<system_t const &>())
      )
     >::value
   and
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
	     std::declval<::pressio::ode::types::step_t const &>(),
	     std::declval<jacobian_t &>()
	    )
     )
   >::value
   >
 > : std::true_type{};
//------------------------------------------------------------------

template<typename T, typename ... args>
using implicit_euler_jacobian_policy =
  implicit_jacobian_policy<
  T, ::pressio::ode::implicitmethods::Euler, 1, args...>;

template<typename T, typename ... args>
using implicit_bdf2_jacobian_policy =
  implicit_jacobian_policy<
  T, ::pressio::ode::implicitmethods::BDF2,2, args...>;
//------------------------------------------------------------------

}}} // namespace pressio::ode::concepts
#endif
