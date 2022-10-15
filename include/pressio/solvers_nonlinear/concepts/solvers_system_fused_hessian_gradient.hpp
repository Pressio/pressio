/*
//@HEADER
// ************************************************************************
//
// solvers_system_fused_hessian_gradient.hpp
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

#ifndef SOLVERS_NONLINEAR_CONSTRAINTS_SOLVERS_SYSTEM_FUSED_HESSIAN_GRADIENT_HPP_
#define SOLVERS_NONLINEAR_CONSTRAINTS_SOLVERS_SYSTEM_FUSED_HESSIAN_GRADIENT_HPP_

namespace pressio{ namespace nonlinearsolvers{

#ifdef PRESSIO_ENABLE_CXX20
template <class T>
concept SystemWithFusedHessianAndGradient =
  /*
    required nested aliases
  */
  requires(){
    typename T::state_type;
    typename T::gradient_type;
    typename T::hessian_type;
    typename T::residual_norm_type;
  }
  /*
    requirements on the nested aliases
  */
  && ::pressio::is_vector_eigen<typename T::state_type>::value
  && ::pressio::is_vector_eigen<typename T::gradient_type>::value
  && ::pressio::is_dense_matrix_eigen<typename T::hessian_type>::value
  && all_have_traits_and_same_scalar<
	typename T::state_type,
	typename T::gradient_type,
	typename T::hessian_type
     >::value
  /*
    compound requirements
  */
  && requires(const T & A,
	      const typename T::state_type & state,
	      typename T::gradient_type & g,
	      typename T::hessian_type & H,
	      pressio::Norm normKind,
	      typename T::residual_norm_type & rNormValue,
	      bool recomputeJacobian)
  {
    { A.createState()       } -> std::same_as<typename T::state_type>;
    { A.createGradient()    } -> std::same_as<typename T::gradient_type>;
    { A.createHessian()     } -> std::same_as<typename T::hessian_type>;

    { A.residualNorm(state, normKind, rNormValue) } -> std::same_as<void>;

    { A.hessianAndGradient(state, H, g, normKind,
			   rNormValue, recomputeJacobian)  } -> std::same_as<void>;
  };
#endif //PRESSIO_ENABLE_CXX20

}} // end namespace pressio::nonlinearsolvers

#if not defined PRESSIO_ENABLE_CXX20
namespace pressio{ namespace nonlinearsolvers{

template<typename T, typename enable = void>
struct SystemWithFusedHessianAndGradient : std::false_type{};

template<typename T>
struct SystemWithFusedHessianAndGradient<
  T,
 ::pressio::mpl::enable_if_t<
      ::pressio::has_state_typedef<T>::value
   && ::pressio::has_hessian_typedef<T>::value
   && ::pressio::has_gradient_typedef<T>::value
   && ::pressio::has_residual_norm_typedef<T>::value
   //
   && ::pressio::is_vector_eigen<typename T::state_type>::value
   && ::pressio::is_vector_eigen<typename T::gradient_type>::value
   && ::pressio::is_dense_matrix_eigen<typename T::hessian_type>::value
   && all_have_traits_and_same_scalar<
	 typename T::state_type,
	 typename T::gradient_type,
	 typename T::hessian_type
      >::value
   //
   && ::pressio::nonlinearsolvers::has_const_create_state_method_return_result<
     T, typename T::state_type>::value
   && ::pressio::nonlinearsolvers::has_const_create_hessian_method_return_result<
     T, typename T::hessian_type>::value
   && ::pressio::nonlinearsolvers::has_const_create_gradient_method_return_result<
     T, typename T::gradient_type>::value
   //
   && ::pressio::nonlinearsolvers::has_const_hessianandgradient_method_accept_state_result_norm_return_void<
	T, typename T::state_type, typename T::hessian_type,
	typename T::gradient_type, typename T::residual_norm_type>::value
   && ::pressio::nonlinearsolvers::has_const_residualnorm_method_accept_state_norm_return_void<
	T, typename T::state_type, typename T::residual_norm_type>::value
   >
 > : std::true_type{};

}} // end namespace pressio::nonlinearsolvers
#endif

#endif  // SOLVERS_NONLINEAR_CONSTRAINTS_SOLVERS_SYSTEM_FUSED_HESSIAN_GRADIENT_HPP_
