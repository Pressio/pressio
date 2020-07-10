
#ifndef rom_continuous_time_system_PRECONDITIONABLE_ROM_HPP_
#define rom_continuous_time_system_PRECONDITIONABLE_ROM_HPP_

namespace pressio{ namespace rom{ namespace concepts {

template<typename T, typename enable = void>
struct continuous_time_system_preconditionable_rom : std::false_type{};

template<typename T>
struct continuous_time_system_preconditionable_rom<
  T,
  mpl::enable_if_t<
    ::pressio::rom::concepts::continuous_time_system<T>::value and 
    ::pressio::rom::predicates::has_const_apply_preconditioner_method_accept_state_time_result_return_void<
            T, typename T::state_type, typename T::scalar_type, typename T::velocity_type >::value and
    ::pressio::rom::predicates::has_const_apply_preconditioner_method_accept_state_time_result_return_void<
            T, typename T::state_type, typename T::scalar_type, typename T::dense_matrix_type >::value 
    >
  > : std::true_type{};


}}} // namespace pressio::rom::concepts
#endif
