
#ifndef rom_steady_system_HPP_
#define rom_steady_system_HPP_

namespace pressio{ namespace rom{ namespace concepts {

template<typename T, typename enable = void>
struct steady_system : std::false_type{};

template<typename T>
struct steady_system<
  T,
  mpl::enable_if_t<
    ::pressio::containers::predicates::has_scalar_typedef<T>::value and
    ::pressio::ode::predicates::has_state_typedef<T>::value and
    ::pressio::ode::predicates::has_residual_typedef<T>::value and
    ::pressio::rom::predicates::has_dense_matrix_typedef<T>::value and
    ///////////////////
    /// residual 
    ///////////////////
    ::pressio::rom::predicates::has_const_create_residual_method_return_result<
      T, typename T::residual_type>::value and
    ::pressio::rom::predicates::has_const_residual_method_accept_state_result_return_void<
      T, typename T::state_type, typename T::residual_type
      >::value and 
    ///////////////////
    /// apply jacobian
    ///////////////////
    ::pressio::rom::predicates::has_const_create_apply_jacobian_result_method_accept_operand_return_result<
      T, typename T::dense_matrix_type, typename T::dense_matrix_type >::value and
    ::pressio::rom::predicates::has_const_apply_jacobian_method_accept_state_operand_result_return_void<
      T, typename T::state_type,  typename T::dense_matrix_type, typename T::dense_matrix_type
      >::value
    >
  > : std::true_type{};


}}} // namespace pressio::rom::concepts
#endif
