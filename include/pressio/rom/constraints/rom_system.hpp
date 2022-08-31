/*
//@HEADER
// ************************************************************************
//
// rom_fom_system_continuous_time.hpp
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef ROM_CONSTRAINTS_ROM_FOM_SYSTEM_HPP_
#define ROM_CONSTRAINTS_ROM_FOM_SYSTEM_HPP_

namespace pressio{ namespace rom{

//
// rhs only
//
template<class T, class enable = void>
struct SemiDiscreteFom : std::false_type{};

template<class T>
struct SemiDiscreteFom<
  T,
  mpl::enable_if_t<
       ::pressio::has_time_typedef<T>::value
    && ::pressio::has_state_typedef<T>::value
    && ::pressio::has_right_hand_side_typedef<T>::value
    //
    && std::is_copy_constructible<typename T::state_type>::value
    && std::is_copy_constructible<typename T::right_hand_side_type>::value
    && ::pressio::VectorSpaceElementsWithSameField<
      typename T::state_type, typename T::right_hand_side_type
      >::value
    && std::is_convertible<
      typename T::time_type,
      typename ::pressio::Traits<typename T::state_type>::scalar_type
      >::value
    //
    && ::pressio::ode::has_const_create_rhs_method_return_result<
      T, typename T::right_hand_side_type >::value
    && ::pressio::ode::has_const_rhs_method_accept_state_indvar_result_return_void<
      T, typename T::state_type, typename T::time_type, typename T::right_hand_side_type>::value
   >
  > : std::true_type{};

//
// rhs, jacobian action
// need to detect what is the type of the jacobian action
//
template<class T, class ManifoldJacType, class enable = void>
struct SemiDiscreteFomWithJacobianAction : std::false_type{};

template<class T, class ManifoldJacType>
struct SemiDiscreteFomWithJacobianAction<
  T, ManifoldJacType,
  mpl::enable_if_t<
       SemiDiscreteFom<T>::value
    //
    && ::pressio::rom::has_const_create_apply_jacobian_result_method_accept_operand_return_result<
	 T, ManifoldJacType>::value
    && ::pressio::rom::has_const_apply_jacobian_method_accept_state_operand_time_result_return_void<
	 T, typename T::state_type, ManifoldJacType, typename T::time_type,
	 // use decltype to deduce the return type of the jac action method
	 decltype(
		  std::declval<T const>().createApplyJacobianResult
		  (
		   std::declval<ManifoldJacType const &>()
		   )
		  )
	 >::value
    && std::is_copy_constructible<
	 decltype(
		  std::declval<T const>().createApplyJacobianResult
		  (
		   std::declval<ManifoldJacType const &>()
		   )
		  )
	 >::value
    && ::pressio::VectorSpaceElementsWithSameField<
	 typename T::state_type,
         decltype(
		  std::declval<T const>().createApplyJacobianResult
		  (
		   std::declval<ManifoldJacType const &>()
		   )
		  )
	 >::value
    >
  > : std::true_type{};


//
// steady fom
//
template<class T, class ManifoldJacType, class enable = void>
struct SteadyFomWithJacobianAction : std::false_type{};

template<class T, class ManifoldJacType>
struct SteadyFomWithJacobianAction<
  T, ManifoldJacType,
  mpl::enable_if_t<
       ::pressio::has_state_typedef<T>::value
    && ::pressio::has_residual_typedef<T>::value
    //
    && std::is_copy_constructible<typename T::state_type>::value
    && std::is_copy_constructible<typename T::residual_type>::value
    && ::pressio::VectorSpaceElementsWithSameField<
      typename T::state_type, typename T::residual_type
      >::value
    //
    && ::pressio::rom::has_const_create_residual_method_return_result<
      T, typename T::residual_type >::value
    && ::pressio::rom::has_const_residual_method_accept_state_result_return_void<
      T, typename T::state_type, typename T::residual_type>::value
    && ::pressio::rom::has_const_create_apply_jacobian_result_method_accept_operand_return_result<
	 T, ManifoldJacType>::value
    && ::pressio::rom::has_const_apply_jacobian_method_accept_state_operand_result_return_void<
	 T, typename T::state_type, ManifoldJacType,
	 // use decltype to deduce the return type of the jac action method
	 decltype(
		  std::declval<T const>().createApplyJacobianResult
		  (
		   std::declval<ManifoldJacType const &>()
		   )
		  )
	 >::value
    && std::is_copy_constructible<
	 decltype(
		  std::declval<T const>().createApplyJacobianResult
		  (
		   std::declval<ManifoldJacType const &>()
		   )
		  )
	 >::value
    && ::pressio::VectorSpaceElementsWithSameField<
	 typename T::state_type,
         decltype(
		  std::declval<T const>().createApplyJacobianResult
		  (
		   std::declval<ManifoldJacType const &>()
		   )
		  )
	 >::value
   >
  > : std::true_type{};

//
// fully discrete
//
template<class T, int NumStates, class ManifoldJacType, class enable = void>
struct FullyDiscreteSystemWithJacobianAction : std::false_type{};

template<class T, int NumStates, class ManifoldJacType>
struct FullyDiscreteSystemWithJacobianAction<
  T, NumStates, ManifoldJacType,
  mpl::enable_if_t<
       ::pressio::has_time_typedef<T>::value
    && ::pressio::has_state_typedef<T>::value
    && ::pressio::has_discrete_residual_typedef<T>::value
    //
    && mpl::is_same<
	 typename T::discrete_residual_type,
	 decltype(std::declval<T const>().createDiscreteTimeResidual())
	 >::value
    //
    // && std::is_void<
    // 	decltype
    // 	(
    // 	 std::declval<T const>().createResultOfDiscreteTimeJacobianAction
    // 	 (
    // 	  std::declval<ManifoldJacType const &>()
    // 	  )
    // 	 )
    // 	>::value
    // //
    // && ::pressio::rom::has_const_discrete_residual_jacobian_action_method<
    // 	 T, NumStates, int,
    // 	 typename T::independent_variable_type,
    // 	 typename T::state_type,
    // 	 typename T::discrete_residual_type,
    // 	 ManifoldJacType,
    // 	 decltype
    // 	 (
    // 	  std::declval<T const>().createResultOfDiscreteTimeJacobianAction
    // 	  (
    // 	   std::declval<ManifoldJacType const &>()
    // 	   )
    // 	 )
    // 	 >::value
    >
  > : std::true_type{};

}}
#endif  // ROM_CONSTRAINTS_ROM_FOM_SYSTEM_CONTINUOUS_TIME_HPP_



// //
// // rhs, mass matrix action
// // need to detect what is the type of the MM action
// //
// template<class T, class ManifoldJacType, class enable = void>
// struct SemiDiscreteFomWithMassMatrixAction : std::false_type{};

// template<class T, class ManifoldJacType>
// struct SemiDiscreteFomWithMassMatrixAction<
//   T, ManifoldJacType,
//   mpl::enable_if_t<
//        SemiDiscreteFom<T>::value
//     //
//     && ::pressio::rom::has_const_create_apply_mass_matrix_result_method_accept_operand_return_result<
// 	 T, ManifoldJacType>::value
//     && ::pressio::rom::has_const_apply_mass_matrix_method_accept_state_operand_time_result_return_void<
// 	 T, typename T::state_type, ManifoldJacType, typename T::time_type,
// 	 // use decltype to deduce the return type of the jac action method
// 	 decltype(
// 		  std::declval<T const>().createApplyMassMatrixResult(
// 								      std::declval<ManifoldJacType const &>()
// 								      )
// 		  )
// 	 >::value
//     && std::is_copy_constructible<
// 	 decltype(
// 		  std::declval<T const>().createApplyMassMatrixResult(
// 								      std::declval<ManifoldJacType const &>()
// 								      )
// 		  )
// 	 >::value
//     && ::pressio::VectorSpaceElementsWithSameField<
// 	 typename T::state_type,
//          decltype(
// 		  std::declval<T const>().createApplyMassMatrixResult(
// 								      std::declval<ManifoldJacType const &>()
// 								      )
// 		  )
// 	 >::value
//     >
//   > : std::true_type{};

// //
// // rhs, constant mass matrix action
// // need to detect what is the type of the MM action
// //
// template<class T, class ManifoldJacType, class enable = void>
// struct SemiDiscreteFomWithConstantMassMatrixAction : std::false_type{};

// template<class T, class ManifoldJacType>
// struct SemiDiscreteFomWithConstantMassMatrixAction<
//   T, ManifoldJacType,
//   mpl::enable_if_t<
//        SemiDiscreteFom<T>::value
//     //
//     && ::pressio::rom::has_const_create_apply_mass_matrix_result_method_accept_operand_return_result<
// 	 T, ManifoldJacType>::value
//     && ::pressio::rom::has_const_apply_mass_matrix_method_accept_operand_result_return_void<
// 	 T, ManifoldJacType,
// 	 // use decltype to deduce the return type of the jac action method
// 	 decltype(
// 		  std::declval<T const>().createApplyMassMatrixResult(
// 								      std::declval<ManifoldJacType const &>()
// 								      )
// 		  )
// 	 >::value
//     && std::is_copy_constructible<
// 	 decltype(
// 		  std::declval<T const>().createApplyMassMatrixResult(
// 								      std::declval<ManifoldJacType const &>()
// 								      )
// 		  )
// 	 >::value
//     && ::pressio::VectorSpaceElementsWithSameField<
// 	 typename T::state_type,
//          decltype(
// 		  std::declval<T const>().createApplyMassMatrixResult(
// 								      std::declval<ManifoldJacType const &>()
// 								      )
// 		  )
// 	 >::value
//     >
//   > : std::true_type{};


// //
// // rhs, mass matrix action
// // need to detect what is the type of the MM action
// //
// template<class T, class ManifoldJacType, class enable = void>
// struct SemiDiscreteFomComplete : std::false_type{};

// template<class T, class ManifoldJacType>
// struct SemiDiscreteFomComplete<
//   T, ManifoldJacType,
//   mpl::enable_if_t<
//        SemiDiscreteFomWithMassMatrixAction<T, ManifoldJacType>::value
//     && SemiDiscreteFomWithJacobianAction<T, ManifoldJacType>::value
//     >
//   > : std::true_type{};
