/*
//@HEADER
// ************************************************************************
//
// rom_galerkin_fom_velocity_policy.hpp
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

#ifndef ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_POLICIES_ROM_GALERKIN_FOM_VELOCITY_POLICY_HPP_
#define ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_POLICIES_ROM_GALERKIN_FOM_VELOCITY_POLICY_HPP_

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

template <typename fom_states_manager_t, typename fom_velocity_type>
class FomVelocityPolicy
{
public:
  using data_type = fom_velocity_type;

private:
  std::reference_wrapper<fom_states_manager_t> fomStatesMngr_;
  mutable fom_velocity_type fomVelo_ = {};

public:
  FomVelocityPolicy() = delete;
  FomVelocityPolicy(const FomVelocityPolicy &) = default;
  FomVelocityPolicy & operator=(const FomVelocityPolicy &) = delete;
  FomVelocityPolicy(FomVelocityPolicy &&) = default;
  FomVelocityPolicy & operator=(FomVelocityPolicy &&) = delete;
  ~FomVelocityPolicy() = default;

  template<typename fom_system_t>
  FomVelocityPolicy(const fom_system_t & fomSystemObj,
		    fom_states_manager_t & fomStatesMngr)
    : fomStatesMngr_(fomStatesMngr),
      fomVelo_(fomSystemObj.createVelocity())
  {}

public:
  template<
  class galerkin_state_t, class fom_system_t, class scalar_t, class at_tag
  >
  void compute(const galerkin_state_t & galerkinState,
	       const fom_system_t  & fomSystemObj,
	       const scalar_t & velocityEvalTime,
	       at_tag tag) const
  {
    fomStatesMngr_.get().reconstructAt(galerkinState, tag);
    const auto & fomState = fomStatesMngr_.get().fomStateAt(tag);
    fomSystemObj.velocity(*fomState.data(), velocityEvalTime, *fomVelo_.data());
  }

  const fom_velocity_type & get() const{
    return fomVelo_;
  }
};

}}}}//end namespace pressio::rom::galerkin::impl
#endif  // ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_POLICIES_ROM_GALERKIN_FOM_VELOCITY_POLICY_HPP_
