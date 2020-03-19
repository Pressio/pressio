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

#ifndef OPTIMIZERS_PARAMS_HPP_
#define OPTIMIZERS_PARAMS_HPP_

namespace pressio{ namespace optimizers{

template <typename scalar_type>
class Parameters{

private:
  scalar_type gradientOptimalityTol_   = {}/* default value */;
  scalar_type stepOptimalityTol_       = {}/* default value */;
  scalar_type constraintOptimalityTol_ = {}/* default value */;
  int maxIters_                      = {}/* default value */;
  int metricsVerbosity_		       = {}/* default value */;

public:
  void setGradientOptimalityTolerance(scalar_type newTol){
    /* minimum objective gradient magnitude at which to terminate */
    gradientOptimalityTol_ = newTol;
  }

  void setObjFunctionChangeOptimalityTolerance(scalar_type newTol){
    /* minimum CHANGE in objective function at which to terminate*/
    stepOptimalityTol_ = newTol;
  }

  void setConstraintOptimalityTolerance(scalar_type newTol){
    //  minimum constraint violation magnitude at which to terminate
    // Conditional: on the gradient magnitude also below its tolerance.
    constraintOptimalityTol_ = newTol;
  }

  void setMaxIterations(int maxIter){
    // maximum number of optimization iterations.
    maxIters_ = maxIter;
  }

  // getters corresponding to above setters
  // -----------

  // Optimization sub problem parameters

  // Diagnostics
  void setInfoVerbosity(int value){
    /* various levels of verbosity
       if value=0: print nothing
       if value=1: print only obj function value
       if value=2: print obj function and gradient value
       etc. */
    metricsVerbosity_ = value;
  }
};

}}//end namespace pressio::optimizers
#endif
