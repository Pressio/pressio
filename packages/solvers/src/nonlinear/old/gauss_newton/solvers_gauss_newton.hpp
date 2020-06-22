/*
//@HEADER
// ************************************************************************
//
// solvers_gauss_newton.hpp
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

#ifndef SOLVERS_GAUSS_NEWTON_HPP_
#define SOLVERS_GAUSS_NEWTON_HPP_

#include "./impl_neq/solvers_gn_neq_specializer_impl.hpp"
#include "./impl_neq/res_jac_api/solvers_gn_neq_py_impl.hpp"
#include "./impl_qr/solvers_gn_qr_specializer_impl.hpp"

namespace pressio{ namespace solvers{ namespace nonlinear{

/* alias: GN solvers normal equations */
template <typename ... Args>
using GaussNewton = typename impl::GaussNewtonNormalEquationsSpecializer<Args...>::type;

/* alias: QR-based GN solvers */
template <typename ... Args>
using GaussNewtonQR = typename impl::GaussNewtonQRSpecializer<Args...>::type;

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
template <typename ... Args>
using PyGaussNewton = typename impl::PyGaussNewton<Args...>::type;
#endif

}}}//end namespace pressio::solvers::nonlinear
#endif
