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
void _printMetrics(bool printGradient,
		   step_t iStep,
		   const sc_t & absCorrectionNorm,
		   const sc_t & relCorrectionNorm,
		   const sc_t & absResNorm,
		   const sc_t & relResNorm,
		   const sc_t & absGNorm,
		   const sc_t & relGNorm)
{
  if (printGradient)
  {
    PRESSIOLOG_INFO
      ("nonlinIter = {:2d}: ||R||(a) = {:.6e} ||R||(r) = {:.6e} ||g||(a) = {:.6e} ||g||(r) = {:.6e} ||delta||(a) = {:.6e} ||delta||(r) = {:.6e}",
       iStep, absResNorm, relResNorm, absGNorm, relGNorm, absCorrectionNorm, relCorrectionNorm);
  }
  else
  {
    PRESSIOLOG_INFO
      ("nonlinIter = {:2d}: ||R||(a) = {:.6e} ||R||(r) = {:.6e} ||delta||(a) = {:.6e} ||delta||(r) = {:.6e}",
       iStep, absResNorm, relResNorm, absCorrectionNorm, relCorrectionNorm);
  }
}

template <typename step_t, typename sc_t>
void printMetrics(step_t iStep,
		  const sc_t & absCorrectionNorm,
		  const sc_t & relCorrectionNorm,
		  const sc_t & absResNorm,
		  const sc_t & relResNorm,
		  const sc_t & absGNorm,
		  const sc_t & relGNorm)
{
  _printMetrics(true,
		iStep, absCorrectionNorm, relCorrectionNorm,
		absResNorm, relResNorm,
		absGNorm, relGNorm);
}

template <typename step_t, typename sc_t>
void printMetrics(step_t iStep,
		  const sc_t & absCorrectionNorm,
		  const sc_t & relCorrectionNorm,
		  const sc_t & absResNorm,
		  const sc_t & relResNorm)
{
  _printMetrics(false,
		iStep, absCorrectionNorm, relCorrectionNorm,
		absResNorm, relResNorm, sc_t(), sc_t());
}

}}}}

#endif  // SOLVERS_NONLINEAR_IMPL_SOLVERS_PRINTER_HPP_
