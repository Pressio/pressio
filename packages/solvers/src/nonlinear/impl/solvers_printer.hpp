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

#ifndef SOLVERS_NONLINEAR_IMPL_PRINTER_SOLVERS_PRINTER_
#define SOLVERS_NONLINEAR_IMPL_PRINTER_SOLVERS_PRINTER_

namespace pressio{ namespace solvers{ namespace nonlinear{ namespace impl{
template <typename sc_t>
class NonlinearLeastSquaresDefaultMetricsPrinter{
private:
  sc_t gNorm0_ = {};
  sc_t resNorm0_ = {};

public:
  template <typename solver_t, typename step_t>
  void givenGradientNormsPrintRest(const solver_t & solver, step_t iStep, const sc_t & absGNorm, const sc_t & relGNorm)
  {   
      const auto correctionNorm = solver.correctionNormCurrentCorrectionStep();
      const auto resNorm = solver.residualNormCurrentCorrectionStep();
      if (iStep == 1) resNorm0_ = resNorm;
      printImpl(iStep, correctionNorm, resNorm, resNorm/resNorm0_, absGNorm, relGNorm);
  }

  template <typename solver_t, typename step_t>
  void givenCorrectionNormsPrintRest(const solver_t & solver, step_t iStep, const sc_t & absCorrectionNorm, const sc_t & relCorrectionNorm)
  {   
      const auto absGNorm = solver.gradientNormCurrentCorrectionStep();
      const auto absResNorm = solver.residualNormCurrentCorrectionStep();
      if (iStep == 1){
        gNorm0_ = absGNorm;
        resNorm0_ = absResNorm;
      }
      const auto relResNorm = absResNorm/resNorm0_;
      const auto relGNorm = absGNorm/gNorm0_;

      printImpl(iStep, absCorrectionNorm, absResNorm, relResNorm, absGNorm, relGNorm);
  }



  template <typename solver_t, typename step_t>
  void givenResidualNormsPrintRest(const solver_t & solver, step_t iStep, const sc_t & absResNorm, const sc_t & relResNorm)
  {   
      const auto correctionNorm = solver.correctionNormCurrentCorrectionStep();
      const auto absGNorm = solver.gradientNormCurrentCorrectionStep();
      if (iStep == 1) gNorm0_ = absGNorm;
      printImpl(iStep, correctionNorm, absResNorm, relResNorm, absGNorm, absGNorm/gNorm0_);
  }

  template <typename solver_t, typename step_t>
  void print(const solver_t & solver, step_t iStep)
  {   
      const auto correctionNorm = solver.correctionNormCurrentCorrectionStep();
      const auto absGNorm = solver.gradientNormCurrentCorrectionStep();
      const auto absResNorm = solver.residualNormCurrentCorrectionStep();
      if (iStep == 1) gNorm0_ = absGNorm;
      if (iStep == 1) resNorm0_ = absResNorm;

      printImpl(iStep, correctionNorm, absResNorm, absResNorm/resNorm0_, absGNorm, absGNorm/gNorm0_);
  }


private:
  template <typename step_t>
  void printImpl(step_t iStep, 
                 const sc_t & correctionNorm,
                 const sc_t & absResNorm, const sc_t & relResNorm, 
                 const sc_t & absGNorm,   const sc_t & relGNorm) const 
  {
      auto fmt = utils::io::cyan() + utils::io::bold();      
      ::pressio::utils::io::print_stdout(fmt, std::scientific,
            " Nonlinear iteration no. =" , iStep,
            utils::io::reset(),
            " Residual l2 norm (abs) =", absResNorm,
            " Residual l2 norm (rel) =", relResNorm,
            " Gradient l2 norm (abs) =", absGNorm,
            " Gradient l2 norm (rel) =", relGNorm,
            " Correction l2 norm =", correctionNorm,
            "\n");
  }

};

}}}}

#endif //SOLVERS_NONLINEAR_IMPL_PRINTER_SOLVERS_PRINTER_

