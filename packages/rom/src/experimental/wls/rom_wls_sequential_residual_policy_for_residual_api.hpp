/*
//@HEADER
// ************************************************************************
//
// rom_wls_sequential_residual_policy_for_residual_api.hpp
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

#ifndef ROM_WLS_SEQUENTIAL_RESIDUAL_POLICY_FOR_RESIDUAL_API_HPP_
#define ROM_WLS_SEQUENTIAL_RESIDUAL_POLICY_FOR_RESIDUAL_API_HPP_

namespace pressio{ namespace rom{ namespace experimental{

template <std::size_t numAuxStates, typename residual_type, typename fom_states_container_type>
class WlsSequentialResidualPolicyForResidualApi
{

public:
  using this_t = WlsSequentialResidualPolicyForResidualApi<residual_type, fom_states_container_type>;
  using residual_t = residual_type;

public:
  WlsSequentialResidualPolicyForResidualApi() = delete;
  ~WlsSequentialResidualPolicyForResidualApi() = default;

  WlsSequentialResidualPolicyForResidualApi(fom_states_container_type & fomStatesIn)
    : fomStates_(fomStatesIn){}

public:
  template <typename wls_state_t, typename fom_t>
  void operator()(const wls_state_t			& wlsState,
  		  const fom_t				& fomObj,
  		  residual_t				& wlsR) const
  {
    this->compute_impl(wlsState, fomObj, wlsR);
  }

  template <typename wls_state_t, typename fom_t>
  residual_t operator()(const wls_state_t		    & wlsState,
			const fom_t			    & fomObj) const
  {
    // here I need to:
    // 1. reconstruct the fom state (this method is typically only called once to init
    //		the data, so it does not really matter which state we reconstruct.
    //		What we care is the shape of the object, not the data in it.
    //		Let's say we pick the fom state at 0. So we reconstruct the fom state at step 0.)
    // 2. query the fom object to get a time-discrete residual object.
    // 3. construct and return R

    // by doing this, we expect the fomState_ to be of a type supporting this method
    // see for example WlsFomStatesContainer.
    fomStates_.template reconstructFomStateAt(wlsState, 0);

    // we need to make sure residual_t (which for wls is probably a container with many instances
    // of the fom time-discrete residual) has a constructor compatible with the following syntax
    // since the fomObj will return a single residual object
    residual_t R( fomObj.createTimeDiscreteResidualObject( *fomStates_[0].data() ));

    return R;
  }


private:
  template <typename wls_state_t, typename fom_t>
  void compute_impl(const wls_state_t		        & wlsState,
		    const fom_t			        & fomObj,
		    residual_t			        & wlsR) const
  {
    /* here is where I need to do things sequentially.
     * at each step in the window, I need to recontruct the fom states at and compute residuals
     * I can do this in two ways:
     * (a) reconstruct all fom states and then compute residuals
     * (b) loop over each step in the window, reconstruct local fom states and compute residual here
     *
     * approach (b) allows us to never have to store all the fom states so we just need
    */

    // just for simplicity for now, let's do (a)
    for (std::size_t iStep=0; iStep < fom_state_container_type.n_; ++iStep){
      fomStates_.reconstructFomStateAt(wlsState, iStep);
    }

    // compute the residual at step 0 in a special way?
    wlsR[0] = ...;

    // now that the fom states have been reconstructed at all steps in the window
    // compute the wls residual by computing the time discrete residual at each step
    for (std::size_t iStep=1; iStep < fom_state_container_type.n_; ++iStep)
    {
      const auto & currentFomState = fomStates_[iStep];
      auto & currentFomResidual = wlsR[iStep];

      // this if should be done at compile time but leave it here for now just to get things down
      if (numAuxStates == 1){
	fomObj.template timeDiscreteResidual(step, time, dt,
					     *currentFomResidual.data(),
					     *currentFomState.data(),
					     *fomStates_[iStep-1].data());
      }
    }
  }

protected:
  fom_states_container_type & fomStates_;
};

}}}//end namespace pressio::rom::experimental
#endif
