
#ifndef ode_continuous_time_api_explicit_ode_HPP_
#define ode_continuous_time_api_explicit_ode_HPP_

namespace pressio{ namespace ode{ namespace concepts {

template<typename T, typename enable = void>
struct continuous_time_system_explicit_stepping : std::false_type{};


template<typename T>
struct continuous_time_system_explicit_stepping<
  T,
  mpl::enable_if_t<
    ::pressio::containers::predicates::has_scalar_typedef<T>::value and
    ::pressio::ode::predicates::has_state_typedef<T>::value and
    ::pressio::ode::predicates::has_velocity_typedef<T>::value and
    //// velocity methods
    ::pressio::ode::predicates::has_const_create_velocity_method_return_result<
      T, typename T::velocity_type >::value and
    ::pressio::ode::predicates::has_const_velocity_method_accept_state_time_result_return_void<
      T, typename T::state_type,typename T::scalar_type,typename T::velocity_type>::value
   >
  > : std::true_type{};


#ifdef PRESSIO_ENABLE_TPL_PYBIND11
template<typename T>
struct continuous_time_system_explicit_stepping<
  T,
  mpl::enable_if_t<
    ::pressio::mpl::is_same<
      T, pybind11::object>::value
    >
  > : std::true_type{};
#endif

}}} // namespace pressio::ode::concepts
#endif
