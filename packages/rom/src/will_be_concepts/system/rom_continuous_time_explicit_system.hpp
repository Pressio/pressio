
#ifndef rom_continuous_time_explicit_system_HPP_
#define rom_continuous_time_explicit_system_HPP_

namespace pressio{ namespace rom{ namespace concepts {

template<typename T, typename enable = void>
struct continuous_time_explicit_system : std::false_type{};

// #ifdef PRESSIO_ENABLE_TPL_PYBIND11
// template<typename T>
// struct continuous_time_explicit_system<
//   T,
//   mpl::enable_if_t<
//     ::pressio::mpl::is_same<T, pybind11::object>::value
//     >
//   > : std::true_type{};
// #endif


template<typename T>
struct continuous_time_explicit_system<
  T,
  mpl::enable_if_t<
    ::pressio::containers::predicates::has_scalar_typedef<T>::value and
    ::pressio::ode::predicates::has_state_typedef<T>::value and
    ::pressio::ode::predicates::has_velocity_typedef<T>::value and
    /// velocity 
    ::pressio::ode::predicates::has_const_create_velocity_method_return_result<
      T, typename T::velocity_type>::value and
    ::pressio::ode::predicates::has_const_velocity_method_accept_state_time_result_return_void<
      T, typename T::state_type, typename T::scalar_type, typename T::velocity_type
      >::value
    >
  > : std::true_type{};

}}} // namespace pressio::rom::concepts
#endif
