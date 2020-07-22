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
template <typename sc_t>
class NonlinearLeastSquaresDefaultMetricsPrinter
{

private:
  sc_t residualNorm0_ = {};
  sc_t gradientNorm0_ = {};

public:
  template <typename solver_t, typename step_t>
  void print(const solver_t & solver, step_t iStep,
	     const sc_t & absoluteCorrecNorm,   const sc_t & relativeCorrecNorm,
	     const sc_t & absoluteResidualNorm,	const sc_t & relativeResidualNorm,
	     const sc_t & absoluteGradientNorm, const sc_t & relativeGradientNorm)
  {
    if (solver.computesGradient()){
      printImpl(iStep,
		absoluteCorrecNorm,
		absoluteResidualNorm, relativeResidualNorm,
		absoluteGradientNorm, relativeGradientNorm);
    }
    else{
      printImpl(iStep, absoluteCorrecNorm,
		absoluteResidualNorm, relativeResidualNorm);
    }
  }

private:
  template <typename step_t>
  void printImpl(step_t iStep,
                 const sc_t & correctionNorm,
                 const sc_t & absResNorm, const sc_t & relResNorm,
                 const sc_t & absGNorm,   const sc_t & relGNorm) const
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
	 "||Correction||_l2 =", correctionNorm,
	 reset(),
	 "\n");
  }

  template <typename step_t>
  void printImpl(step_t iStep, const sc_t & correctionNorm,
                 const sc_t & absResNorm, const sc_t & relResNorm) const
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
	 "||Correction||_l2 =", correctionNorm,
	 reset(),
	 "\n");
  }

  // template <typename solver_t, typename step_t>
  // void print(const solver_t & solver, step_t iStep)
  // {
  //   const auto correctionNorm	   = solver.correctionNormCurrentCorrectionStep();
  //   const auto absResNorm	   = solver.residualNormCurrentCorrectionStep();
  //   if (iStep == 1) residualNorm0_ = absResNorm;

  //   if (solver.computesGradient()){
  //     const auto absGNorm	      = solver.gradientNormCurrentCorrectionStep();
  //     if (iStep == 1) gradientNorm0_  = absGNorm;
  //     printImpl(iStep, correctionNorm, absResNorm, absResNorm/residualNorm0_,
  // 		absGNorm, absGNorm/gradientNorm0_);
  //   }
  //   else{
  //     printImpl(iStep, correctionNorm, absResNorm, absResNorm/residualNorm0_);
  //   }
  // }

  // template <typename solver_t, typename step_t>
  // void givenGradientNormsPrintRest(const solver_t & solver,
  // 				   step_t iStep,
  // 				   const sc_t & absGNorm,
  // 				   const sc_t & relGNorm)
  // {
  //   const auto correctionNorm = solver.correctionNormCurrentCorrectionStep();
  //   const auto resNorm	      = solver.residualNormCurrentCorrectionStep();

  //   if (iStep == 1) residualNorm0_ = resNorm;
  //   printImpl(iStep, correctionNorm, resNorm, resNorm/residualNorm0_,
  // 	      absGNorm, relGNorm);
  // }

  // template <typename solver_t, typename step_t>
  // void givenCorrectionNormsPrintRest(const solver_t & solver,
  // 				     step_t iStep,
  // 				     const sc_t & absCorrectionNorm,
  // 				     const sc_t & relCorrectionNorm)
  // {
  //   const auto absResNorm = solver.residualNormCurrentCorrectionStep();
  //   if (iStep == 1) residualNorm0_ = absResNorm;

  //   if (solver.computesGradient()){
  //     const auto absGNorm   = solver.gradientNormCurrentCorrectionStep();
  //     if (iStep == 1) gradientNorm0_  = absGNorm;
  //     printImpl(iStep, absCorrectionNorm, absResNorm, absResNorm/residualNorm0_,
  // 		absGNorm, absGNorm/gradientNorm0_);
  //   }
  //   else{
  //     printImpl(iStep, absCorrectionNorm, absResNorm, absResNorm/residualNorm0_);
  //   }
  // }

  // template <typename solver_t, typename step_t>
  // void givenResidualNormsPrintRest(const solver_t & solver,
  // 				   step_t iStep,
  // 				   const sc_t & absResNorm,
  // 				   const sc_t & relResNorm)
  // {
  //   const auto correctionNorm = solver.correctionNormCurrentCorrectionStep();

  //   if (solver.computesGradient()){
  //     const auto absGNorm	      = solver.gradientNormCurrentCorrectionStep();
  //     if (iStep == 1) gradientNorm0_  = absGNorm;
  //     printImpl(iStep, correctionNorm, absResNorm, relResNorm,
  // 		absGNorm, absGNorm/gradientNorm0_);
  //   }
  //   else{
  //     printImpl(iStep, correctionNorm, absResNorm, relResNorm);
  //   }
  // }

};

}}}}

#endif  // SOLVERS_NONLINEAR_IMPL_SOLVERS_PRINTER_HPP_
