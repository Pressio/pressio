
#ifndef rom_admissible_system_time_discrete_residual_api_unsteady_lspg_hpp_
#define rom_admissible_system_time_discrete_residual_api_unsteady_lspg_hpp_

namespace pressio{ namespace rom{ namespace meta {

template<typename T, typename enable = void>
struct admissible_system_time_discrete_residual_api_unsteady_lspg : std::false_type{};

// #ifdef PRESSIO_ENABLE_TPL_PYBIND11
// template<typename T>
// struct rom_admissible_system_time_discrete_residual_api_unsteady_lspg<
//   T,
//   mpl::enable_if_t<
//     ::pressio::mpl::is_same<T, pybind11::object>::value
//     >
//   > : std::true_type{};
// #endif

template<typename T>
struct admissible_system_time_discrete_residual_api_unsteady_lspg<
  T,
  mpl::enable_if_t<
    ::pressio::containers::predicates::has_scalar_typedef<T>::value and
    ::pressio::ode::predicates::has_state_typedef<T>::value and
    ::pressio::ode::predicates::has_time_discrete_residual_typedef<T>::value and
    ::pressio::rom::meta::has_dense_matrix_typedef<T>::value 
    and
    ///////////////////////////
    // time-discrete residual 
    ::pressio::ode::predicates::has_const_create_td_residual_method_return_result<
        T, typename T::time_discrete_residual_type>::value 
    and 
    ::pressio::ode::predicates::has_const_td_residual_method_accept_step_time_dt_result_norm_n_states_return_void<
        T, 1, ::pressio::ode::types::step_t, 
        typename T::scalar_type, 
        typename T::state_type, 
        typename T::time_discrete_residual_type>::value 
    and 
    ::pressio::ode::predicates::has_const_td_residual_method_accept_step_time_dt_result_norm_n_states_return_void<
        T, 2, ::pressio::ode::types::step_t, 
        typename T::scalar_type, 
        typename T::state_type, 
        typename T::time_discrete_residual_type>::value 
    and
    ::pressio::ode::predicates::has_const_td_residual_method_accept_step_time_dt_result_norm_n_states_return_void<
        T, 3, ::pressio::ode::types::step_t, 
        typename T::scalar_type, 
        typename T::state_type, 
        typename T::time_discrete_residual_type>::value
    and        
    // 
    ///////////////////////////
    // apply time-discrete jacobian
    has_const_create_apply_td_jacobian_result_method_accept_operand_return_result<
        T,  typename T::dense_matrix_type, typename T::dense_matrix_type>::value 
    and 
    has_const_apply_td_jacobian_method_accept_step_time_dt_operand_result_n_states_returning_void<
        T, 1, 
        ::pressio::ode::types::step_t, 
        typename T::scalar_type, 
        typename T::state_type, 
        typename T::dense_matrix_type,
        typename T::dense_matrix_type>::value 
    and 
    has_const_apply_td_jacobian_method_accept_step_time_dt_operand_result_n_states_returning_void<
        T, 2, 
        ::pressio::ode::types::step_t, 
        typename T::scalar_type, 
        typename T::state_type, 
        typename T::dense_matrix_type,
        typename T::dense_matrix_type>::value 
    and 
    has_const_apply_td_jacobian_method_accept_step_time_dt_operand_result_n_states_returning_void<
        T, 3, 
        ::pressio::ode::types::step_t, 
        typename T::scalar_type, 
        typename T::state_type, 
        typename T::dense_matrix_type,
        typename T::dense_matrix_type>::value
    >
  > : std::true_type{};


}}} // namespace pressio::rom::meta
#endif
