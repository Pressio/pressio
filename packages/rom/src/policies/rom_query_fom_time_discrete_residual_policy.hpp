/*
//@HEADER
// ************************************************************************
//
// rom_query_fom_time_discrete_residual_policy.hpp
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

#ifndef ROM_QUERY_FOM_TIME_DISCRETE_RESIDUAL_HPP_
#define ROM_QUERY_FOM_TIME_DISCRETE_RESIDUAL_HPP_

namespace pressio{ namespace rom{ namespace policy{

struct QueryFomTimeDiscreteResidual
{

  // ---------------------
  // for n = 1
  // ---------------------
  template <typename fom_t, typename step_t, typename time_t, typename result_t, typename fom_state_t>
  void evaluate(const fom_state_t & fomCurrentState,
  		const std::array<fom_state_t, 1> & fomPrevStates,
		const fom_t   & fomObj,
		const time_t  & time,
		const time_t  & dt,
  		const step_t  & step,
  		result_t      & R) const
  {
    fomObj.template timeDiscreteResidual(step, time, dt, *R.data(),
    					 *fomCurrentState.data(),
    					 *fomPrevStates[0].data());
  }

  template <typename fom_t, typename step_t, typename time_t, typename fom_state_t>
  auto evaluate(const fom_state_t & fomCurrentState,
		const std::array<fom_state_t, 1> & fomPrevStates,
		const fom_t   & fomObj,
  		const time_t  & time,
  		const time_t  & dt,
  		const step_t  & step) const
    -> decltype(
  		fomObj.template timeDiscreteResidual(step, time, dt,
						     *fomCurrentState.data(),
						     *(fomPrevStates[0].data()))
  		)
  {
    return fomObj.template timeDiscreteResidual(step, time, dt, *fomCurrentState.data(),
						*(fomPrevStates[0].data()));
  }


  // ---------------------
  // for n = 2
  // ---------------------
  template <typename fom_t, typename step_t, typename time_t, typename result_t, typename fom_state_t>
  void evaluate(const fom_state_t & fomCurrentState,
  		const std::array<fom_state_t, 2> & fomPrevStates,
		const fom_t   & fomObj,
		const time_t  & time,
		const time_t  & dt,
  		const step_t  & step,
  		result_t      & R) const
  {
    fomObj.template timeDiscreteResidual(step, time, dt, *R.data(),
    					 *fomCurrentState.data(),
    					 *fomPrevStates[0].data(),
					 *fomPrevStates[1].data());
  }

  template <typename fom_t, typename step_t, typename time_t, typename fom_state_t>
  auto evaluate(const fom_state_t & fomCurrentState,
		const std::array<fom_state_t, 2> & fomPrevStates,
		const fom_t   & fomObj,
  		const time_t  & time,
  		const time_t  & dt,
  		const step_t  & step) const
    -> decltype(
  		fomObj.template timeDiscreteResidual(step, time, dt,
						     *fomCurrentState.data(),
						     *fomPrevStates[0].data(),
						     *fomPrevStates[1].data())
  		)
  {
    return fomObj.template timeDiscreteResidual(step, time, dt,
						*fomCurrentState.data(),
						*fomPrevStates[0].data(),
						*fomPrevStates[1].data());
  }

};

}}} //end namespace pressio::rom::policy
#endif
