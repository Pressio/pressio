

#ifndef ODE_ode_admissible_system_implicit_ode_arbitrary_stepper_HPP_
#define ODE_ode_admissible_system_implicit_ode_arbitrary_stepper_HPP_

namespace pressio{ namespace ode{ namespace meta {

template<typename T, typename enable = void>
struct admissible_system_implicit_ode_arbitrary_stepper : std::false_type{};

template<typename T>
struct admissible_system_implicit_ode_arbitrary_stepper<
  T,
  mpl::enable_if_t<
    ::pressio::containers::meta::has_scalar_typedef<T>::value and
    ::pressio::ode::meta::has_state_typedef<T>::value and
    ::pressio::ode::meta::has_residual_typedef<T>::value and
    ::pressio::ode::meta::has_jacobian_typedef<T>::value and
    //
    // time-discrete residual 
    has_const_create_td_residual_method_return_result<
        T, typename T::residual_type>::value and 
    has_const_td_residual_method_accept_step_time_dt_result_norm_states_return_void<
        T, 1, types::step_t, typename T::scalar_type, typename T::scalar_type, 
        typename T::residual_type>::value and 
    has_const_td_residual_method_accept_step_time_dt_result_norm_states_return_void<
        T, 2, types::step_t, typename T::scalar_type, typename T::scalar_type, 
        typename T::residual_type>::value and
    has_const_td_residual_method_accept_step_time_dt_result_norm_states_return_void<
        T, 3, types::step_t, typename T::scalar_type, typename T::scalar_type, 
        typename T::residual_type>::value and        
    // 
    // time-discrete jacobian
    has_const_create_td_jacobian_method_return_result<
        T, typename T::jacobian_type>::value and 
    has_const_td_jacobian_method_accept_step_time_dt_result_states_return_void<
        T, 1, types::step_t, typename T::scalar_type, typename T::scalar_type, 
        typename T::jacobian_type>::value and 
    has_const_td_jacobian_method_accept_step_time_dt_result_states_return_void<
        T, 2, types::step_t, typename T::scalar_type, typename T::scalar_type, 
        typename T::jacobian_type>::value and 
    has_const_td_jacobian_method_accept_step_time_dt_result_states_return_void<
        T, 3, types::step_t, typename T::scalar_type, typename T::scalar_type, 
        typename T::jacobian_type>::value 
    >
  > : std::true_type{};


#ifdef PRESSIO_ENABLE_TPL_PYBIND11
template<typename T>
struct admissible_system_implicit_ode_arbitrary_stepper<
  T,
  mpl::enable_if_t<
    mpl::is_same<T, pybind11::object>::value
    >
  > : std::true_type{};
#endif

}}} // namespace pressio::ode::meta
#endif
