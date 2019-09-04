/*
//@HEADER
// ************************************************************************
//
// solvers_convergence_tags.hpp
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

#ifndef SOLVERS_CONVERGENCE_TAGS_HPP_
#define SOLVERS_CONVERGENCE_TAGS_HPP_

#include "solvers_norm_tags.hpp"

namespace pressio{ namespace solvers{ namespace iterative{

namespace converged_when{

template <typename norm_type>
struct absoluteNormCorrectionBelowTol{
  using norm_t = norm_type;
  static_assert(mpl::is_same<norm_t, ::pressio::solvers::L1Norm>::value or
		mpl::is_same<norm_t, ::pressio::solvers::L2Norm>::value,
		"Invalid template for absoluteNormCorrectionBelowTol, it must be L1Norm or L2Norm");
};

template <typename norm_type>
struct absoluteNormResidualBelowTol{
  using norm_t = norm_type;
  static_assert(mpl::is_same<norm_t, ::pressio::solvers::L1Norm>::value or
		mpl::is_same<norm_t, ::pressio::solvers::L2Norm>::value,
		"Invalid template for absoluteNormResidualBelowTol, it must be L1Norm or L2Norm");
};

template <typename norm_type>
struct relativeNormResidualBelowTol{
  using norm_t = norm_type;
  static_assert(mpl::is_same<norm_t, ::pressio::solvers::L1Norm>::value or
		mpl::is_same<norm_t, ::pressio::solvers::L2Norm>::value,
		"Invalid template for relativeNormResidualBelowTol, it must be L1Norm or L2Norm");
};

template <typename norm_type>
struct absoluteNormProjectedResidualBelowTol{
  using norm_t = norm_type;
  static_assert(mpl::is_same<norm_t, ::pressio::solvers::L1Norm>::value or
		mpl::is_same<norm_t, ::pressio::solvers::L2Norm>::value,
		"Invalid template for absoluteNormProjectedResidualBelowTol, it must be L1Norm or L2Norm");
};

template <typename norm_type>
struct relativeNormProjectedResidualBelowTol{
  using norm_t = norm_type;
  static_assert(mpl::is_same<norm_t, ::pressio::solvers::L1Norm>::value or
		mpl::is_same<norm_t, ::pressio::solvers::L2Norm>::value,
		"Invalid template for relativeNormProjectedResidualBelowTol, it must be L1Norm or L2Norm");
};


struct completingNumMaxIters{};

}//end namespace pressio::solvers::convergedWhen

using default_convergence
	= converged_when::absoluteNormCorrectionBelowTol<L2Norm>;

}}}//end namespace pressio::solvers::iterative
#endif
