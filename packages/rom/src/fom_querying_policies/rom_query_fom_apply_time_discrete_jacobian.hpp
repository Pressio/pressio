/*
//@HEADER
// ************************************************************************
//
// rom_query_fom_apply_time_discrete_jacobian_policy.hpp
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

#ifndef ROM_QUERY_FOM_APPLY_TIME_DISCRETE_JACOBIAN_HPP_
#define ROM_QUERY_FOM_APPLY_TIME_DISCRETE_JACOBIAN_HPP_

namespace pressio{ namespace rom{ namespace policy{

struct QueryFomApplyTimeDiscreteJacobian
{

  template <class fom_state_t, class fom_t, class operand_t>
  auto evaluate(const fom_state_t & fomCurrentState,
  		const fom_t	  & fomObj,
  		const operand_t   & B) const
    -> decltype(
  		fomObj.createApplyTimeDiscreteJacobianObject(*fomCurrentState.data(),
							     *B.data())
  		)
  {
    return fomObj.createApplyTimeDiscreteJacobianObject(*fomCurrentState.data(),
							*B.data());
  }

  template <
    class fom_state_t, class fom_t, class step_t, class time_t, class operand_t, class result_t
    >
  void evaluate(const fom_state_t & state_n,
		const fom_state_t & state_nm1,
  		const fom_t	  & fomObj,
  		const time_t	  & time,
  		const time_t	  & dt,
  		const step_t	  & step,
  		const operand_t   & B,
  		result_t	  & A,
		// by default, we compute wrt current state
  		int compute_jac_wrt_state_id = 0) const
  {
    fomObj.template applyTimeDiscreteJacobian(step, time, dt,
  					      *B.data(),
  					      compute_jac_wrt_state_id,
  					      *A.data(),
					      *state_n.data(),
					      *state_nm1.data());
  }


  template <
    class fom_state_t, class fom_t, class step_t, class time_t, class operand_t, class result_t
    >
  void evaluate(const fom_state_t & state_n,
		const fom_state_t & state_nm1,
		const fom_state_t & state_nm2,
  		const fom_t	  & fomObj,
  		const time_t	  & time,
  		const time_t	  & dt,
  		const step_t	  & step,
  		const operand_t   & B,
  		result_t	  & A,
		// by default, we compute wrt current state
  		int compute_jac_wrt_state_id = 0) const
  {
    fomObj.template applyTimeDiscreteJacobian(step, time, dt,
  					      *B.data(),
  					      compute_jac_wrt_state_id,
  					      *A.data(),
					      *state_n.data(),
					      *state_nm1.data(),
					      *state_nm2.data());
  }
};

}}} //end namespace pressio::rom::policy
#endif
