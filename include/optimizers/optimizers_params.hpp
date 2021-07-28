/*
//@HEADER
// ************************************************************************
//
// optimizers_params.hpp
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

#ifndef OPTIMIZERS_OPTIMIZERS_PARAMS_HPP_
#define OPTIMIZERS_OPTIMIZERS_PARAMS_HPP_

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "ROL_ParameterList.hpp"
#endif

namespace pressio{ namespace optimizers{

template <typename scalar_type>
class Parameters{

private:
  ::pressio::optimizers::stepMethod chosenStepMethod_ = default_step_method;
  int maxIters_			       = 100;
  scalar_type gradientOptimalityTol_   = 1e-10;
  scalar_type stepOptimalityTol_       = 1e-12;
  //  int metricsVerbosity_	       = {}/* default value */;

public:
  void setStepMethod(::pressio::optimizers::stepMethod stepMethodIn){
    chosenStepMethod_ = stepMethodIn;
  }
  ::pressio::optimizers::stepMethod getStepMethod() const{
    return chosenStepMethod_;
  }

  void setMaxIterations(int maxIter){
    // maximum number of optimization iterations.
    maxIters_ = maxIter;
  }
  int maxIterations() const{
    return maxIters_;
  }

  void setGradientNormOptimalityTolerance(scalar_type newTol){
    /* minimum objective gradient magnitude at which to terminate */
    gradientOptimalityTol_ = newTol;
  }
  scalar_type getGradientNormOptimalityTolerance() const{
    return gradientOptimalityTol_;
  }

  void setStepNormOptimalityTolerance(scalar_type newTol){
    /* norm of the step at which to terminate*/
    stepOptimalityTol_ = newTol;
  }
  scalar_type getStepNormOptimalityTolerance() const{
    return stepOptimalityTol_;
  }

  // // Diagnostics
  // void setInfoVerbosity(int value){
  //   /* various levels of verbosity
  //      if value=0: print nothing
  //      if value=1: print only obj function value
  //      if value=2: print obj function and gradient value
  //      etc. */
  //   metricsVerbosity_ = value;
  // }
};

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
template <typename scalar_type>
void convertToRolParameterList(const Parameters<scalar_type> & params,
			       ROL::ParameterList & rolParList)
{
  rolParList.sublist("Status Test").set("Iteration Limit", params.maxIterations());
  rolParList.sublist("Status Test").set("Gradient Tolerance", params.getGradientNormOptimalityTolerance());
  rolParList.sublist("Status Test").set("Step Tolerance", params.getStepNormOptimalityTolerance());

  const auto e = params.getStepMethod();
  switch (e)
    {
    case ::pressio::optimizers::stepMethod::lineSearch:
      {
	rolParList.sublist("Step").set("Type","Line Search");
	rolParList.sublist("Step").sublist("Line Search").set("Subproblem Solver","Truncated CG");
	break;
      }
    case ::pressio::optimizers::stepMethod::trustRegion:
      {
	rolParList.sublist("Step").set("Type","Trust Region");
	rolParList.sublist("Step").sublist("Trust Region").set("Subproblem Solver","Truncated CG");
	break;
      }
    default:
      {
	break;
      }
    };
}
#endif

}}//end namespace pressio::optimizers
#endif  // OPTIMIZERS_OPTIMIZERS_PARAMS_HPP_
