
#ifndef rom_continuous_time_system_hpp_
#define rom_continuous_time_system_hpp_

namespace pressio{ namespace rom{ namespace concepts {

template<typename T, typename enable = void>
struct continuous_time_system : std::false_type{};

// #ifdef PRESSIO_ENABLE_TPL_PYBIND11
// template<typename T>
// struct continuous_time_system<
//   T,
//   mpl::enable_if_t<
//     ::pressio::mpl::is_same<T, pybind11::object>::value
//     >
//   > : std::true_type{};
// #endif


template<typename T>
struct continuous_time_system<
  T,
  mpl::enable_if_t<
    ::pressio::containers::predicates::has_scalar_typedef<T>::value and
    ::pressio::ode::predicates::has_state_typedef<T>::value and
    ::pressio::ode::predicates::has_velocity_typedef<T>::value and
    // this is the case when we don't have jacobian and dense matrix
    // we need also these nagative ones otherwise template is ambiguous
    !::pressio::ode::predicates::has_jacobian_typedef<T>::value and
    !::pressio::rom::predicates::has_dense_matrix_typedef<T>::value and
    /// velocity 
    ::pressio::ode::predicates::has_const_create_velocity_method_return_result<
      T, typename T::velocity_type>::value and
    ::pressio::ode::predicates::has_const_velocity_method_accept_state_time_result_return_void<
      T, typename T::state_type, typename T::scalar_type, typename T::velocity_type
      >::value
    >
  > : std::true_type{};


template<typename T>
struct continuous_time_system<
  T,
  mpl::enable_if_t<
    ::pressio::containers::predicates::has_scalar_typedef<T>::value and
    ::pressio::ode::predicates::has_state_typedef<T>::value and
    ::pressio::ode::predicates::has_velocity_typedef<T>::value and
    ::pressio::ode::predicates::has_jacobian_typedef<T>::value and
    ::pressio::rom::predicates::has_dense_matrix_typedef<T>::value and
    ///////////////////
    /// velocity 
    ///////////////////
    ::pressio::ode::predicates::has_const_create_velocity_method_return_result<
      T, typename T::velocity_type>::value and
    ::pressio::ode::predicates::has_const_velocity_method_accept_state_time_result_return_void<
      T, typename T::state_type, typename T::scalar_type, typename T::velocity_type
      >::value and 
    ///////////////////
    /// apply jacobian
    ///////////////////
    ::pressio::rom::predicates::has_const_create_apply_jacobian_result_method_accept_operand_return_result<
      T, typename T::dense_matrix_type, typename T::dense_matrix_type >::value and
    ::pressio::rom::predicates::has_const_apply_jacobian_method_accept_state_operand_time_result_return_void<
      T, 
      typename T::state_type,  typename T::dense_matrix_type,
      typename T::scalar_type, typename T::dense_matrix_type
      >::value
    >
  > : std::true_type{};


}}} // namespace pressio::rom::concepts
#endif
