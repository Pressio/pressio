/*
//@HEADER
// ************************************************************************
//
// solvers_fwd.hpp
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

#ifndef SOLVERS_FORWARD_DECLARATIONS_HPP_
#define SOLVERS_FORWARD_DECLARATIONS_HPP_

namespace pressio{ namespace solvers{

namespace direct{
template<typename SolverT, typename MatrixT, typename enable = void>
class KokkosDirect;
}//end namespace pressio::solvers::direct

namespace linear { namespace details {
template <typename T>
struct traits;
}}//end namespace pressio::solvers::linear::details


namespace iterative{ namespace impl{

namespace experimental{
template <
  typename system_type,
  typename linear_solver_type,
  typename scalar_type,
  typename line_search_type,
  typename convergence_when_t
  >
class GaussNewtonHessianGradientApi;
}//end namespace experimental

template <
  typename system_t,
  typename hessian_t,
  typename linear_solver_t,
  typename scalar_t,
  typename line_search_t,
  typename when_converged_t,
  typename resid_obs_t,
  typename ud_ops_t,
  typename enable = void
  >
class GaussNewtonNormalEqResJacApi;

template <
  typename system_t,
  typename qr_solver_t,
  typename scalar_t,
  typename line_search_t,
  typename when_converged_t,
  typename enable = void
  >
class GaussNewtonQR;

template <typename ... Args>
struct GNNEQSpecializationPicker;

template <typename ... Args>
struct GNQRSpecializationPicker;

}//end namespace pressio::solvers::iterative::impl


/* alias: GN solvers normal equations */
template <typename ... Args>
using GaussNewton = typename impl::GNNEQSpecializationPicker<Args...>::type;


/* alias: QR-based GN solvers */
template <typename ... Args>
using GaussNewtonQR = typename impl::GNQRSpecializationPicker<Args...>::type;

/* class to interface with python */
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
template <
  typename system_t,
  typename state_t,
  typename residual_t,
  typename jacobian_t,
  typename hessian_t,
  typename linear_solver_t,
  typename scalar_t,
  typename ops_t,
  typename when_converged_t = ::pressio::solvers::iterative::default_convergence,
  typename enable = void
  >
struct PyGaussNewton;
#endif

namespace hacked{
template <
  typename scalar_t,
  typename lin_solver_tag,
  template <typename, typename> class lin_solver_t,
  typename line_search_t,
  typename when_converged_t = ::pressio::solvers::iterative::default_convergence,
  typename system_t = void,
  typename cbar_t = void,
  typename enable = void
  >
class GaussNewtonConservative;
}//end namespace pressio::solvers::iterative::hacked

}//end namespace pressio::solvers::iterative

}}//end namespace pressio::solvers

#endif
