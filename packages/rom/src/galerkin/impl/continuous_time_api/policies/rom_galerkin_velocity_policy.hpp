/*
//@HEADER
// ************************************************************************
//
// rom_galerkin_explicit_velocity_policy.hpp
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

#ifndef ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_POLICIES_ROM_GALERKIN_VELOCITY_POLICY_HPP_
#define ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_POLICIES_ROM_GALERKIN_VELOCITY_POLICY_HPP_

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

template <typename galerkin_velocity_t, typename projection_policy_t>
class VelocityPolicy : private projection_policy_t
{
  static_assert
  (::pressio::rom::galerkin::concepts::galerkin_velocity<galerkin_velocity_t>::value,
   "The galerkin_velocity_t is not a valid galerkin velocity type");

public:
  VelocityPolicy() = delete;
  VelocityPolicy(const VelocityPolicy &) = default;
  VelocityPolicy & operator=(const VelocityPolicy &) = delete;
  VelocityPolicy(VelocityPolicy &&) = default;
  VelocityPolicy & operator=(VelocityPolicy &&) = delete;
  ~VelocityPolicy() = default;

  template<typename ...Args>
  VelocityPolicy(std::size_t romSize,
		 Args && ...args)
    : projection_policy_t(std::forward<Args>(args)...),
      romSize_(romSize){}

public:
  template <typename fom_system_t>
  galerkin_velocity_t create(const fom_system_t & fomSystemObj) const
  {
    galerkin_velocity_t result(romSize_);
    ::pressio::ops::set_zero(result);
    return result;
  }

  template<class galerkin_state_t, typename fom_system_t, typename scalar_t>
  void compute(const galerkin_state_t & galerkinState,
	       galerkin_velocity_t & galerkinRhs,
	       const fom_system_t  & fomSystemObj,
	       const scalar_t & time) const
  {
    projection_policy_t::compute(galerkinRhs, galerkinState, fomSystemObj, time);
  }

protected:
  const std::size_t romSize_ = {};
};

}}}}//end namespace pressio::rom::galerkin::impl
#endif  // ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_POLICIES_ROM_GALERKIN_VELOCITY_POLICY_HPP_
