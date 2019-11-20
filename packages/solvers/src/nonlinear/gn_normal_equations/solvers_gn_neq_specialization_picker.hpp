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

#include "../../solvers_fwd.hpp"
#include "../../meta/solvers_is_legitimate_system_for_nonlinear_solver.hpp"
#include "../../meta/solvers_is_legitimate_linear_solver_for_gn_normeq.hpp"
#include "../../meta/solvers_is_legitimate_hessian_for_gn_normeq.hpp"
#include "../../meta/solvers_is_legitimate_line_search_tag.hpp"
#include "../../meta/solvers_is_legitimate_convergence_tag.hpp"
#include "../../meta/solvers_is_legitimate_residual_observer_when_converged.hpp"

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
  typename system_t, typename scalar_t, typename linear_solver_t,
  typename line_search_t, typename convergence_t, typename ... Args>
struct GNNEQWithResidualJacobainApi
{
  using linear_solver_matrix_t = typename linear_solver_t::matrix_type;

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

  // check if sequence contains an observer for the residual after GN is converged
  using ic6 = ::pressio::mpl::variadic::find_if_binary_pred_t<
    typename system_t::residual_type,
    ::pressio::solvers::meta::is_legitimate_residual_observer_when_solver_converged, Args...>;
  // store the type
  using observer_when_conv_t = ::pressio::mpl::variadic::at_or_t<void, ic6::value, Args...>;

  // check if sequence contains an observer for the residual at each GN step
  using ic7 = ::pressio::mpl::variadic::find_if_binary_pred_t<
    typename system_t::residual_type,
    ::pressio::solvers::meta::is_legitimate_residual_observer_each_solver_step, Args...>;
  // store the type
  using observer_each_step_t = ::pressio::mpl::variadic::at_or_t<void, ic7::value, Args...>;

  /* currently we only allows all observation methods on the residual to be
   * in the same observer class. So we need to check here that when both
   *  observers type detected are NOT void, then they must be the same type */
  using obs_supported_t = ObserverTypesSupported<observer_when_conv_t, observer_each_step_t>;
  static_assert( obs_supported_t::value,
  "Currently, can only accept a single residual observer type to do everything. \
  You can still have multiple methods inside that type to cover all possible observations,\
  but cannot pass multiple types with individual observation methods");
  using observer_t = typename obs_supported_t::type;

  using type = ::pressio::solvers::iterative::impl::GaussNewton<
    system_t, hessian_t, linear_solver_t, scalar_t, line_search_t, convergence_t, observer_t>;
};



template <
  typename system_t, typename scalar_t, typename linear_solver_t,
  typename line_search_t, typename convergence_t, typename ... Args>
struct GNNEQWithHessianGradientApi
{
  using type = ::pressio::solvers::iterative::impl::experimental::GaussNewtonHessianGradientApi<
    system_t, linear_solver_t, scalar_t, line_search_t, convergence_t>;
};




template <typename ... Args>
struct GNNEQSpecializationPicker{

  /* ------------------------------------------------ */
  // verify the sequence contains a valid system type
  using ic1 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::solvers::meta::is_legitimate_system_for_nonlinear_solver, Args...>;
  using system_t = ::pressio::mpl::variadic::at_or_t<void, ic1::value, Args...>;
  static_assert(!std::is_void<system_t>::value and ic1::value < sizeof... (Args),
		"A valid system type must be passed to GN templates");

  /* ------------------------------------------------ */
  // since system is valid, detect the scalar type
  using scalar_t = typename system_t::scalar_type;

  /* ------------------------------------------------ */
  // verify the sequence contains a valid solver type
  using ic2 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::solvers::meta::is_legitimate_linear_solver_for_gn_normeq, Args...>;
  using linear_solver_t = ::pressio::mpl::variadic::at_or_t<void, ic2::value, Args...>;
  static_assert(!std::is_void<linear_solver_t>::value and ic2::value < sizeof... (Args),
  		"A valid linear solver type must be passed to GN templates");

  /* ------------------------------------------------ */
  // check if sequence contains a line search tag
  using ic4 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::solvers::meta::is_legitimate_line_search_tag, Args...>;
  using default_no_ls = ::pressio::solvers::iterative::gn::noLineSearch;
  using line_search_t = ::pressio::mpl::variadic::at_or_t<default_no_ls, ic4::value, Args...>;
  static_assert(!std::is_void<line_search_t>::value,
		"The line search type for GN cannot be void");

  /* ------------------------------------------------ */
  // check if sequence contains a valid default convergence
  using ic5 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::solvers::meta::is_legitimate_convergence_tag, Args...>;
  using default_conv = ::pressio::solvers::iterative::default_convergence;
  using convergence_t = ::pressio::mpl::variadic::at_or_t<default_conv, ic5::value, Args...>;
  static_assert(!std::is_void<convergence_t>::value, "The convergence type for GN cannot be void");

  /* ------------------------------------------------ */
  // the types above are common for all APIs, not pass args...
  // to specializers for further inspection
  using type = typename std::conditional<
    ::pressio::solvers::meta::experimental::is_legitimate_system_for_gn_hessian_gradient_api<system_t>::value,
    typename GNNEQWithHessianGradientApi<system_t, scalar_t, linear_solver_t, line_search_t, convergence_t, Args...>::type,
    typename GNNEQWithResidualJacobainApi<system_t, scalar_t, linear_solver_t, line_search_t, convergence_t, Args...>::type
    >::type;
};

}}}}//end namespace pressio::solvers::iterative::impl
#endif
