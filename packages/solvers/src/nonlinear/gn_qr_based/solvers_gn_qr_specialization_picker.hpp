/*
//@HEADER
// ************************************************************************
//
// solvers_gn_qr_specialization_picker.hpp
//                     		      Pressio 
// Copyright 2019 National Technology & Engineering Solutions of Sandia,LLC 
//							      (NTESS)
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

#ifndef SOLVERS_GN_QR_SPECIALIZATION_PICKER_HPP_
#define SOLVERS_GN_QR_SPECIALIZATION_PICKER_HPP_

#include "../../solvers_fwd.hpp"
#include "../../meta/solvers_is_legitimate_system_for_nonlinear_solver.hpp"
#include "../../meta/solvers_is_legitimate_qr_solver_for_gn_qr.hpp"
#include "../../meta/solvers_is_legitimate_line_search_tag.hpp"
#include "../../meta/solvers_is_legitimate_convergence_tag.hpp"

namespace pressio{ namespace solvers{ namespace iterative{ namespace impl{

template <typename ... Args>
struct GNQRSpecializationPicker{

  // verify the sequence contains a valid system type
  using ic1 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::solvers::meta::is_legitimate_system_for_nonlinear_solver, Args...>;
  // store the type
  using system_t = ::pressio::mpl::variadic::at_or_t<void, ic1::value, Args...>;
  static_assert(!std::is_void<system_t>::value and
		ic1::value < sizeof... (Args),
		"A valid system type must be passed to GN templates");


  // verify the sequence contains a valid QR solver type
  using ic2 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::solvers::meta::is_legitimate_qr_solver_for_gn_qr, Args...>;
  // store the type
  using qr_solver_t = ::pressio::mpl::variadic::at_or_t<void, ic2::value, Args...>;
  static_assert(!std::is_void<qr_solver_t>::value and
		ic2::value < sizeof... (Args),
  		"A valid QR solver type must be passed to GN templates");
  // get the type of the matrix from the QR solver
  using qr_solver_matrix_t = typename ::pressio::qr::details::traits<qr_solver_t>::matrix_t;


  // check if sequence contains a legitimate line search tag
  using ic4 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::solvers::meta::is_legitimate_line_search_tag, Args...>;
  // store the type
  using default_no_ls = ::pressio::solvers::iterative::gn::noLineSearch;
  using line_search_t = ::pressio::mpl::variadic::at_or_t<default_no_ls, ic4::value, Args...>;
  static_assert(!std::is_void<line_search_t>::value,
  		"The line search type for GN cannot be void");


  // check if sequence contains a valid default convergence
  using ic5 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::solvers::meta::is_legitimate_convergence_tag, Args...>;
  // store the type
  using default_conv = ::pressio::solvers::iterative::default_convergence;
  using convergence_t = ::pressio::mpl::variadic::at_or_t<default_conv, ic5::value, Args...>;
  static_assert(!std::is_void<convergence_t>::value,
  		"The convergence type for GN cannot be void");

  // using types detected above, define type of GN solver implementation
  using scalar_t = typename system_t::scalar_type;

  using type = ::pressio::solvers::iterative::impl::GaussNewtonQR<
    system_t, qr_solver_t, scalar_t, line_search_t, convergence_t>;
};

}}}}//end namespace pressio::solvers::iterative::impl
#endif
