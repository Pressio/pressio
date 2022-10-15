/*
//@HEADER
// ************************************************************************
//
// solvers_system_residual_jacobian.hpp
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

#ifndef SOLVERS_NONLINEAR_CONSTRAINTS_SOLVERS_SYSTEM_RESIDUAL_JACOBIAN_HPP_
#define SOLVERS_NONLINEAR_CONSTRAINTS_SOLVERS_SYSTEM_RESIDUAL_JACOBIAN_HPP_

namespace pressio{ namespace nonlinearsolvers{

#ifdef PRESSIO_ENABLE_CXX20
template <class T>
concept SystemWithResidualAndJacobian =
  /*
    required nested aliases
  */
  requires(){
    typename T::state_type;
    typename T::residual_type;
    typename T::jacobian_type;
  }
  /*
    requirements on the nested aliases
  */
  && ::pressio::ops::is_known_data_type<typename T::state_type>::value
  && ::pressio::ops::is_known_data_type<typename T::residual_type>::value
  && ::pressio::ops::is_known_data_type<typename T::jacobian_type>::value
  && all_have_traits_and_same_scalar<
      typename T::state_type, typename T::residual_type, typename T::jacobian_type
     >::value
  /*
    compound requirements
  */
  && requires(const T & A,
	      const typename T::state_type & state,
	      typename T::residual_type & r,
	      typename T::jacobian_type & j)
  {
    { A.createState()       } -> std::same_as<typename T::state_type>;
    { A.createResidual()    } -> std::same_as<typename T::residual_type>;
    { A.createJacobian()    } -> std::same_as<typename T::jacobian_type>;
    { A.residual(state, r)  } -> std::same_as<void>;
    { A.jacobian(state, j)  } -> std::same_as<void>;
  };
#endif //PRESSIO_ENABLE_CXX20

}} // end namespace pressio::nonlinearsolvers


namespace pressio{ namespace nonlinearsolvers{

#ifdef PRESSIO_ENABLE_CXX20
template <class T>
concept DeterminedSystemWithResidualAndJacobian = SystemWithResidualAndJacobian<T>;
#endif //PRESSIO_ENABLE_CXX20

}} // end namespace pressio::nonlinearsolvers


namespace pressio{ namespace nonlinearsolvers{

#ifdef PRESSIO_ENABLE_CXX20
template <class T>
concept OverdeterminedSystemWithResidualAndJacobian = SystemWithResidualAndJacobian<T>;
#endif //PRESSIO_ENABLE_CXX20

}} // end namespace pressio::nonlinearsolvers


#if not defined PRESSIO_ENABLE_CXX20

namespace pressio{ namespace nonlinearsolvers{

template<class T, class enable = void>
struct SystemWithResidualAndJacobian : std::false_type{};

template<class T>
struct SystemWithResidualAndJacobian<
  T,
  mpl::enable_if_t<
       ::pressio::has_state_typedef<T>::value
    && ::pressio::has_residual_typedef<T>::value
    && ::pressio::has_jacobian_typedef<T>::value
    //
    && ::pressio::ops::is_known_data_type<typename T::state_type>::value
    && ::pressio::ops::is_known_data_type<typename T::residual_type>::value
    && ::pressio::ops::is_known_data_type<typename T::jacobian_type>::value
    && all_have_traits_and_same_scalar<
	 typename T::state_type,
	 typename T::residual_type,
	 typename T::jacobian_type>::value
    //
    && ::pressio::nonlinearsolvers::has_const_create_state_method_return_result<
      T, typename T::state_type>::value
    && ::pressio::nonlinearsolvers::has_const_create_residual_method_return_result<
      T, typename T::residual_type>::value
    && ::pressio::nonlinearsolvers::has_const_create_jacobian_method_return_result<
      T, typename T::jacobian_type>::value
    //
    && ::pressio::nonlinearsolvers::has_const_residual_method_accept_state_result_return_void<
      T, typename T::state_type, typename T::residual_type>::value
    && ::pressio::nonlinearsolvers::has_const_jacobian_method_accept_state_result_return_void<
      T, typename T::state_type, typename T::jacobian_type>::value
   >
  > : std::true_type{};

template<typename T>
using DeterminedSystemWithResidualAndJacobian = SystemWithResidualAndJacobian<T>;

template<typename T>
using OverdeterminedSystemWithResidualAndJacobian = SystemWithResidualAndJacobian<T>;

}} // end namespace pressio::nonlinearsolvers
#endif // #if not defined PRESSIO_ENABLE_CXX20

#endif  // SOLVERS_NONLINEAR_CONSTRAINTS_SOLVERS_SYSTEM_RESIDUAL_JACOBIAN_HPP_
