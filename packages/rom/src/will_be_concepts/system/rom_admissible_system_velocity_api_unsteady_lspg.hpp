
#ifndef rom_admissible_system_velocity_api_unsteady_lspg_hpp_
#define rom_admissible_system_velocity_api_unsteady_lspg_hpp_

namespace pressio{ namespace rom{ namespace meta {

template<typename T, typename enable = void>
struct admissible_system_velocity_api_unsteady_lspg : std::false_type{};

// #ifdef PRESSIO_ENABLE_TPL_PYBIND11
// template<typename T>
// struct admissible_system_velocity_api_unsteady_lspg<
//   T,
//   mpl::enable_if_t<
//     ::pressio::mpl::is_same<T, pybind11::object>::value
//     >
//   > : std::true_type{};
// #endif

template<typename T>
struct admissible_system_velocity_api_unsteady_lspg<
  T,
  mpl::enable_if_t<
    ::pressio::containers::meta::has_scalar_typedef<T>::value and
    ::pressio::ode::meta::has_state_typedef<T>::value and
    ::pressio::ode::meta::has_velocity_typedef<T>::value and
    ::pressio::ode::meta::has_jacobian_typedef<T>::value and
    ::pressio::rom::meta::has_dense_matrix_typedef<T>::value and
    ///////////////////
    /// velocity 
    ///////////////////
    ::pressio::ode::meta::has_const_create_velocity_method_return_result<
      T, typename T::velocity_type>::value and
    ::pressio::ode::meta::has_const_velocity_method_accept_state_time_result_return_void<
      T, typename T::state_type, typename T::scalar_type, typename T::velocity_type
      >::value and 
    ///////////////////
    /// apply jacobian
    ///////////////////
    has_const_create_apply_jacobian_result_method_accept_operand_return_result<
      T, typename T::dense_matrix_type, typename T::dense_matrix_type >::value and
    has_const_apply_jacobian_method_accept_state_operand_time_result_return_void<
      T, 
      typename T::state_type,  typename T::dense_matrix_type,
      typename T::scalar_type, typename T::dense_matrix_type
      >::value
    >
  > : std::true_type{};


}}} // namespace pressio::rom::meta
#endif
