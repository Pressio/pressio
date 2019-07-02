
#ifndef ODE_IS_LEGITIMATE_AUXILIARY_STEPPER_HPP_
#define ODE_IS_LEGITIMATE_AUXILIARY_STEPPER_HPP_

#include <type_traits>

namespace pressio{ namespace ode{ namespace meta {

template<typename T,
	 typename aux_t,
	 typename enable = void>
struct is_legitimate_auxiliary_stepper : std::false_type{};

template<typename T, typename aux_t>
struct is_legitimate_auxiliary_stepper<
  T, aux_t,
    mpl::enable_if_t<
      mpl::is_same<
	typename ode::details::traits<T>::state_t,
	typename ode::details::traits<aux_t>::state_t
	>::value and
      ::pressio::mpl::is_same<
	typename ode::details::traits<T>::residual_t,
	typename ode::details::traits<aux_t>::residual_t
	>::value and
      ::pressio::mpl::is_same<
	typename ode::details::traits<T>::jacobian_t,
	typename ode::details::traits<aux_t>::jacobian_t
	>::value and
      ::pressio::mpl::is_same<
	typename ode::details::traits<T>::scalar_t,
	typename ode::details::traits<aux_t>::scalar_t
	>::value
      >
  > : std::true_type{};

}}} // namespace pressio::ode::meta
#endif
