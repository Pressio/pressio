
#ifndef rom_continuous_time_system_maskable_rom_HPP_
#define rom_continuous_time_system_maskable_rom_HPP_

namespace pressio{ namespace rom{ namespace concepts {

template<typename T, typename enable = void>
struct continuous_time_system_maskable_rom : std::false_type{};

template<typename T>
struct continuous_time_system_maskable_rom<
  T,
  mpl::enable_if_t<
    ::pressio::rom::concepts::continuous_time_system<T>::value and 
  	// createApplyMaskResult for residual   
    ::pressio::rom::predicates::has_const_create_apply_mask_result_method_accept_operand_return_result<
        T, typename T::velocity_type, typename T::velocity_type>::value and
  	// createApplyMaskResult for dense operator
    ::pressio::rom::predicates::has_const_create_apply_mask_result_method_accept_operand_return_result<
        T, typename T::dense_matrix_type, typename T::dense_matrix_type >::value and
   // applyMask for residual 
    ::pressio::rom::predicates::has_const_apply_mask_method_accept_operand_time_result_return_void<
        T, typename T::velocity_type, typename T::scalar_type, typename T::velocity_type >::value and            
  	// applyMask for dense operator
    ::pressio::rom::predicates::has_const_apply_mask_method_accept_operand_time_result_return_void<
        T, typename T::dense_matrix_type, typename T::scalar_type, typename T::dense_matrix_type >::value 
    >
  > : std::true_type{};


}}} // namespace pressio::rom::concepts
#endif
