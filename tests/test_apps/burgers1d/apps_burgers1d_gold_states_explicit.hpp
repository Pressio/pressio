/*
//@HEADER
// ************************************************************************
//
// apps_burgers1d_gold_states_explicit.hpp
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

#ifndef APPS_BURGERS1D_APPS_BURGERS1D_GOLD_STATES_EXPLICIT_HPP_
#define APPS_BURGERS1D_APPS_BURGERS1D_GOLD_STATES_EXPLICIT_HPP_

namespace pressio { namespace apps{ namespace test{

struct Burgers1dExpGoldStatesEuler{
  using result_t = std::vector<double>;

  static result_t get(int N,
		      double dt,
		      double final_t = 35.,
		      double mu1 = 5.0,
		      double mu2 = 0.02,
		      double mu3 = 0.02)
  {
    if (N == 20 and dt == 0.01 and final_t == 35. and
        mu1 == 5.0 and mu2 == 0.02 and mu3 == 0.02)
    {
      return { 5.0209814000128,  5.044067908724,
	       5.0694601439541,  5.0973757621592,
	       5.1280505161248,  5.1617393082963,
	       5.1987172243105,  5.2392805237326,
	       5.2837475435221,  5.3324594086071,
	       5.3857802812742,  5.4440964817745,
	       5.5078129073313,  5.5773432783592,
	       5.6530870659136,  5.7353794504736,
	       5.8243903774842,  5.9199350492773,
	       6.0211454752168,  6.1259551255163};
    }
    else if (N == 50 and dt == 0.01 and final_t == 35. and
        mu1 == 5.0 and mu2 == 0.02 and mu3 == 0.02)
    {
      return {5.0081549603823, 5.0166286518862, 5.0254329867556,
	      5.0345802812429, 5.0440832662962, 5.0539550983142,
	      5.064209369952, 5.0748601209605, 5.0859218490417,
	      5.0974095206992, 5.1093385820636, 5.1217249696722,
	      5.1345851211777, 5.147935985965, 5.1617950356493,
	      5.1761802744303, 5.1911102492774, 5.2066040599176,
	      5.2226813685989, 5.2393624096019, 5.2566679984713,
	      5.2746195409383, 5.2932390415076, 5.3125491116795,
	      5.3325729777811, 5.3533344883799, 5.3748581212549,
	      5.3971689898997, 5.4202928495376, 5.4442561026257,
	      5.4690858038302, 5.4948096644577, 5.5214560563264,
	      5.5490540150678, 5.577633242846, 5.6072241104857,
	      5.6378576589872, 5.6695656003899, 5.7023803178868,
	      5.7363348649687, 5.7714629630915, 5.8077989967455,
	      5.8453780035106, 5.8842356539985, 5.9244082111713,
	      5.9659324478199, 6.0088454803145, 6.053184437725,
	      6.0989858135867, 6.146284218731};
    }
    else
      return {};

  }//end get
};//end struct

}}}//end namespace pressio::apps::test
#endif  // APPS_BURGERS1D_APPS_BURGERS1D_GOLD_STATES_EXPLICIT_HPP_
