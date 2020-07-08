

#ifndef ode_discrete_time_system_implicit_stepping_HPP_
#define ode_discrete_time_system_implicit_stepping_HPP_

namespace pressio{ namespace ode{ namespace concepts {

template<typename T, typename enable = void>
struct discrete_time_system_implicit_stepping : std::false_type{};

template<typename T>
struct discrete_time_system_implicit_stepping<
  T,
  mpl::enable_if_t<
    ::pressio::containers::predicates::has_scalar_typedef<T>::value and
    ::pressio::ode::predicates::has_state_typedef<T>::value and
    ::pressio::ode::predicates::has_discrete_time_residual_typedef<T>::value and
    ::pressio::ode::predicates::has_discrete_time_jacobian_typedef<T>::value and
    //
    // time-discrete residual 
    ::pressio::ode::predicates::has_const_create_discrete_time_residual_method_return_result<
        T, typename T::discrete_time_residual_type>::value and 
    ::pressio::ode::predicates::has_const_discrete_time_residual_method_accept_step_time_dt_result_norm_n_states_return_void<
        T, 1, types::step_t, 
        typename T::scalar_type, 
        typename T::state_type, 
        typename T::discrete_time_residual_type>::value and 
    ::pressio::ode::predicates::has_const_discrete_time_residual_method_accept_step_time_dt_result_norm_n_states_return_void<
        T, 2, types::step_t, 
        typename T::scalar_type, 
        typename T::state_type, 
        typename T::discrete_time_residual_type>::value and
    // 
    // time-discrete jacobian
    ::pressio::ode::predicates::has_const_create_discrete_time_jacobian_method_return_result<
        T, typename T::discrete_time_jacobian_type>::value and 
    ::pressio::ode::predicates::has_const_discrete_time_jacobian_method_accept_step_time_dt_result_n_states_return_void<
        T, 1, types::step_t, 
        typename T::scalar_type, 
        typename T::state_type, 
        typename T::discrete_time_jacobian_type>::value and 
    ::pressio::ode::predicates::has_const_discrete_time_jacobian_method_accept_step_time_dt_result_n_states_return_void<
        T, 2, types::step_t, 
        typename T::scalar_type, 
        typename T::state_type, 
        typename T::discrete_time_jacobian_type>::value 
    >
  > : std::true_type{};


#ifdef PRESSIO_ENABLE_TPL_PYBIND11
template<typename T>
struct discrete_time_system_implicit_stepping<
  T,
  mpl::enable_if_t<
    mpl::is_same<T, pybind11::object>::value
    >
  > : std::true_type{};
#endif

}}} // namespace pressio::ode::concepts
#endif
