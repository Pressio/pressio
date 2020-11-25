/*
//@HEADER
// ************************************************************************
//
// solvers_printer.hpp
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

#ifndef SOLVERS_NONLINEAR_IMPL_SOLVERS_PRINTER_HPP_
#define SOLVERS_NONLINEAR_IMPL_SOLVERS_PRINTER_HPP_

namespace pressio{ namespace solvers{ namespace nonlinear{ namespace impl{

template <typename step_t, typename sc_t>
void printMetrics(step_t iStep,
		  const sc_t & absCorrectionNorm,
		  const sc_t & relCorrectionNorm,
		  const sc_t & absResNorm,
		  const sc_t & relResNorm,
		  const sc_t & absGNorm,
		  const sc_t & relGNorm)
{
  using namespace ::pressio::utils::io;
  constexpr auto one = static_cast<sc_t>(1);

  // generic format for metrics
  const auto fmt = utils::io::cyan();// + utils::io::bold();

  // use color to highlight if relative residual norm is
  // decreasing (green) or increasing (yellow)
  const auto fmtResRelNorm = relResNorm <= one ?
    (relResNorm < one ? green() : fmt ) : yellow();

  // use color to highlight if relative gradient norm is
  // decreasing (green) or increasing (yellow)
  const auto fmtGradRelNorm = relGNorm <= one ?
    (relGNorm < one ? green() : fmt ) : yellow();

  ::pressio::utils::io::print_stdout
      (
       std::scientific,
       fmt,
       "NonlinearIter =" , iStep,
       " ||Residual||_l2 (abs) =", absResNorm,
       " ||Residual||_l2 (rel) =", fmtResRelNorm, relResNorm, reset(),
       fmt,
       "||Gradient||_l2 (abs) = ", absGNorm,
       " ||Gradient||_l2 (rel) =", fmtGradRelNorm, relGNorm, reset(),
       fmt,
       "||Correction||_l2 =", absCorrectionNorm,
       reset(),
       "\n");
}

template <typename step_t, typename sc_t>
void printMetrics(step_t iStep,
		  const sc_t & absCorrectionNorm,
		  const sc_t & relCorrectionNorm,
		  const sc_t & absResNorm,
		  const sc_t & relResNorm)
{
  using namespace ::pressio::utils::io;
  constexpr auto one = static_cast<sc_t>(1);
  // generic format for metrics
  const auto fmt = utils::io::cyan();// + utils::io::bold();
  // use color to highlight if relative residual norm is
  // decreasing (green) or increasing (yellow)
  const auto fmtResRelNorm = relResNorm <= one ?
    (relResNorm < one ? green() : fmt ) : yellow();

  ::pressio::utils::io::print_stdout
      (std::scientific,
       fmt,
       "NonlinearIter =" , iStep,
       " ||Residual||_l2 (abs) =", absResNorm,
       " ||Residual||_l2 (rel) =", fmtResRelNorm, relResNorm, reset(),
       fmt,
       "||Correction||_l2 =", absCorrectionNorm,
       reset(),
       "\n");
}

}}}}

#endif  // SOLVERS_NONLINEAR_IMPL_SOLVERS_PRINTER_HPP_
