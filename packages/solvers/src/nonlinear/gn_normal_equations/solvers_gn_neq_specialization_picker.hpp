/*
//@HEADER
// ************************************************************************
//
// solvers_gn_neq_specialization_picker.hpp
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

#ifndef SOLVERS_GN_NEQ_SPECIALIZATION_PICKER_HPP_
#define SOLVERS_GN_NEQ_SPECIALIZATION_PICKER_HPP_

namespace pressio{ namespace solvers{ namespace iterative{ namespace impl{

template <typename T1, typename T2>
struct ObserverTypesSupported{
  static constexpr bool value = false;
};

template <>
struct ObserverTypesSupported<void, void>{
  static constexpr bool value = true;
  using type = void;
};

template <typename T>
struct ObserverTypesSupported<void, T>{
  static constexpr bool value = true;
  using type = T;
};

template <typename T>
struct ObserverTypesSupported<T, void>{
  static constexpr bool value = true;
  using type = T;
};

template <typename T>
struct ObserverTypesSupported<T, T>{
  static constexpr bool value = true;
  using type = T;
};




template <
  bool hasDefaultApi,
  typename system_t,
  typename scalar_t,
  typename linear_solver_t,
  typename line_search_t,
  typename convergence_t,
  typename ... Args>
struct GNNEQSpecializeApi;


template <
  typename system_t,
  typename scalar_t,
  typename linear_solver_t,
  typename line_search_t,
  typename convergence_t,
  typename ... Args>
struct GNNEQSpecializeApi<
  true, system_t, scalar_t, linear_solver_t, line_search_t, convergence_t, Args...>
{
  using linear_solver_matrix_t = typename linear_solver_t::matrix_type;

  /* ------------------------------------------------ */
  // check if the sequence contains a valid hessian type
  using ic3 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::solvers::meta::is_legitimate_hessian_for_gn_normeq, Args...>;
  // if no hessian is found, then get it from the linear solver type
  using hessian_t = ::pressio::mpl::variadic::at_or_t<linear_solver_matrix_t, ic3::value, Args...>;
  // in every scenario, the hesssian type must match the matrix type of linear solver
  static_assert( std::is_same<hessian_t, linear_solver_matrix_t>::value,
		 "Hessian type passed to GN must match the matrix type in the linear solver");
  static_assert(!std::is_void<hessian_t>::value,
  		"The hessian type cannot be void");

  /* ------------------------------------------------ */
  // check if sequence contains an observer for the residual after GN is converged
  using ic6 = ::pressio::mpl::variadic::find_if_binary_pred_t<
    typename system_t::residual_type,
    ::pressio::solvers::meta::is_legitimate_residual_observer_when_solver_converged, Args...>;
  using observer_when_conv_t = ::pressio::mpl::variadic::at_or_t<void, ic6::value, Args...>;

  /* ------------------------------------------------ */
  // check if sequence contains an observer for the residual at each GN step
  using ic7 = ::pressio::mpl::variadic::find_if_binary_pred_t<
    typename system_t::residual_type,
    ::pressio::solvers::meta::is_legitimate_residual_observer_each_solver_step, Args...>;
  using observer_each_step_t = ::pressio::mpl::variadic::at_or_t<void, ic7::value, Args...>;

  /* currently we only allow all methods to observe the residual to be
   * in the same observer class. So we need to check here that when both
   *  observers type detected are NOT void, then they must be the same type */
  using obs_supported_t = ObserverTypesSupported<observer_when_conv_t, observer_each_step_t>;
  static_assert( obs_supported_t::value,
  "Currently, can only accept a single residual observer type to do everything. \
  You can still have multiple methods inside that type to cover all possible observations,\
  but cannot pass multiple types with individual observation methods");
  using observer_t = typename obs_supported_t::type;


  /* --------------------------------------------------------- */
  // check if the sequence contains a valid type for custom ops
  /* There are three functionalities needed for GN normal-eq:
   * 1. computing the norm of the residual and gradient
   * 2. computing the hessian J^T * J
   * 3. computing the gradient J^T * R
   *
   * For 2., remember that the system's jacobian can be either a matrix or a mv wrapper
   * and, depending on that, we (might) do things differently when we have to compute hessian.
   * so when we detect custom ops, we need to differentiate based on jacobian's type
   */
  using gradient_t = typename system_t::state_type; // the gradient is same type as state
  using residual_t = typename system_t::residual_type;
  using jacobian_t = typename system_t::jacobian_type;
  using ud_ops_t   = typename GaussNewtonNEqCustomOpsDetectionHelper<gradient_t, residual_t, jacobian_t,
								     hessian_t,  Args...>::type;

  // the class type
  using type = ::pressio::solvers::iterative::impl::GaussNewtonNormalEqResJacApi<
    system_t, hessian_t, linear_solver_t, scalar_t, line_search_t,
    convergence_t, observer_t, ud_ops_t>;
};


template <
  typename system_t,
  typename scalar_t,
  typename linear_solver_t,
  typename line_search_t,
  typename convergence_t,
  typename ... Args>
struct GNNEQSpecializeApi<
  false, system_t, scalar_t, linear_solver_t, line_search_t, convergence_t, Args...>
{
  using type = ::pressio::solvers::iterative::impl::experimental::GaussNewtonHessianGradientApi<
    system_t, linear_solver_t, scalar_t, line_search_t, convergence_t>;
};



template <typename ... Args>
struct GNNEQSpecializationPicker{

  /* ------------------------------------------------ */
  // verify the sequence contains a valid system type
  using ic1 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::solvers::meta::is_legitimate_system_for_gauss_newton_normal_eq, Args...>;
  using system_t = ::pressio::mpl::variadic::at_or_t<void, ic1::value, Args...>;
  static_assert(!std::is_void<system_t>::value and ic1::value < sizeof... (Args),
		"A valid system type must be passed to GN templates. \
This compile-time error means that template arguments passed to the GaussNewton solver\
do not contain a type that is admissible for the normal-equation solver.");
  static constexpr bool hasDefaultApi = ::pressio::solvers::meta::system_meets_default_api<system_t>::value;

  /* ------------------------------------------------ */
  // since system is valid, detect the scalar type
  using scalar_t = typename system_t::scalar_type;

  /* ------------------------------------------------ */
  // check for a valid linear solver type
  using ic2 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::solvers::meta::is_legitimate_linear_solver_for_gn_normeq, Args...>;
  using linear_solver_t = ::pressio::mpl::variadic::at_or_t<void, ic2::value, Args...>;
  static_assert(!std::is_void<linear_solver_t>::value and ic2::value < sizeof... (Args),
  		"A valid linear solver type must be passed to GN with normal equations");

  /* ------------------------------------------------ */
  // check if sequence contains a line search tag
  using ic4 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::solvers::meta::is_legitimate_line_search_tag, Args...>;
  using default_no_ls = ::pressio::solvers::iterative::gn::noLineSearch;
  using line_search_t = ::pressio::mpl::variadic::at_or_t<default_no_ls, ic4::value, Args...>;
  static_assert(!std::is_void<line_search_t>::value,
		"The line search type for GN cannot be void: either omit it so that I use the \
default, or pick one that is valid.");

  /* ------------------------------------------------ */
  // check if sequence contains a valid convergence tag
  using ic5 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::solvers::meta::is_legitimate_convergence_tag, Args...>;
  using default_conv = ::pressio::solvers::iterative::default_convergence;
  using convergence_t = ::pressio::mpl::variadic::at_or_t<default_conv, ic5::value, Args...>;
  static_assert(!std::is_void<convergence_t>::value, "The convergence type for GN cannot be void");

  /* ------------------------------------------------ */
  // the types above are common for all APIs, now pass args to specializers for further inspection
  using type =
    typename GNNEQSpecializeApi<hasDefaultApi, system_t, scalar_t, linear_solver_t, line_search_t, convergence_t, Args...>::type;
};

}}}}//end namespace pressio::solvers::iterative::impl
#endif
