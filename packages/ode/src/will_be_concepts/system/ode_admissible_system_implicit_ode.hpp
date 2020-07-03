
#ifndef ode_admissible_system_implicit_ode_HPP_
#define ode_admissible_system_implicit_ode_HPP_

namespace pressio{ namespace ode{ namespace meta {

template<typename T, typename enable = void>
struct admissible_system_implicit_ode : std::false_type{};

template<typename T>
struct admissible_system_implicit_ode<
  T,
  mpl::enable_if_t<
    admissible_system_implicit_ode_regular_stepper<T>::value
    or
    admissible_system_implicit_ode_arbitrary_stepper<T>::value
    >
  > : std::true_type{};

}}} // namespace pressio::ode::meta
#endif
