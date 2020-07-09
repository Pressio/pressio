
#ifndef rom_discrete_time_system_HPP_
#define rom_discrete_time_system_HPP_

namespace pressio{ namespace rom{ namespace concepts {

template<typename T, typename enable = void>
struct discrete_time_system : std::false_type{};

// #ifdef PRESSIO_ENABLE_TPL_PYBIND11
// template<typename T>
// struct rom_discrete_time_system<
//   T,
//   mpl::enable_if_t<
//     ::pressio::mpl::is_same<T, pybind11::object>::value
//     >
//   > : std::true_type{};
// #endif

template<typename T>
struct discrete_time_system<
  T,
  mpl::enable_if_t<
    ::pressio::containers::predicates::has_scalar_typedef<T>::value and
    ::pressio::ode::predicates::has_state_typedef<T>::value and
    ::pressio::ode::predicates::has_discrete_time_residual_typedef<T>::value and
    ::pressio::rom::predicates::has_dense_matrix_typedef<T>::value 
    and
    ///////////////////////////
    // time-discrete residual 
    ::pressio::ode::predicates::has_const_create_discrete_time_residual_method_return_result<
        T, typename T::discrete_time_residual_type>::value 
    and 
    ::pressio::ode::predicates::has_const_discrete_time_residual_method_accept_step_time_dt_result_norm_n_states_return_void<
        T, 1, ::pressio::ode::types::step_t, 
        typename T::scalar_type, 
        typename T::state_type, 
        typename T::discrete_time_residual_type>::value 
    and 
    ::pressio::ode::predicates::has_const_discrete_time_residual_method_accept_step_time_dt_result_norm_n_states_return_void<
        T, 2, ::pressio::ode::types::step_t, 
        typename T::scalar_type, 
        typename T::state_type, 
        typename T::discrete_time_residual_type>::value 
    and
    // 
    ///////////////////////////
    // apply time-discrete jacobian
    ::pressio::rom::predicates::has_const_create_apply_discrete_time_jacobian_result_method_accept_operand_return_result<
        T,  typename T::dense_matrix_type, typename T::dense_matrix_type>::value 
    and 
    ::pressio::rom::predicates::has_const_apply_discrete_time_jacobian_method_accept_step_time_dt_operand_result_n_states_returning_void<
        T, 1, 
        ::pressio::ode::types::step_t, 
        typename T::scalar_type, 
        typename T::state_type, 
        typename T::dense_matrix_type,
        typename T::dense_matrix_type>::value 
    and 
    ::pressio::rom::predicates::has_const_apply_discrete_time_jacobian_method_accept_step_time_dt_operand_result_n_states_returning_void<
        T, 2, 
        ::pressio::ode::types::step_t, 
        typename T::scalar_type, 
        typename T::state_type, 
        typename T::dense_matrix_type,
        typename T::dense_matrix_type>::value 
    >
  > : std::true_type{};


}}} // namespace pressio::rom::concepts
#endif
