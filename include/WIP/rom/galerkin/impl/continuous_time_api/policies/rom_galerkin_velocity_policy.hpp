/*
//@HEADER
// ************************************************************************
//
// rom_galerkin_velocity_policy.hpp
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
  (::pressio::rom::galerkin::constraints::galerkin_velocity<galerkin_velocity_t>::value,
   "The galerkin_velocity_t is not a valid galerkin velocity type");

  using velocity_traits = ::pressio::containers::details::traits<galerkin_velocity_t>;
  using size_type = typename velocity_traits::size_t;

  // galerkin explicit we can have a rank-1 or rank-2 velocity
  // so when we initialize the velocity here we need to know the rank
  static constexpr int rank_ = velocity_traits::rank;

  // the extents type tshould be changed to something else maybe?
  const std::array<size_type, rank_> extents_ = {};

public:
  VelocityPolicy() = delete;
  VelocityPolicy(const VelocityPolicy &) = default;
  VelocityPolicy & operator=(const VelocityPolicy &) = delete;
  VelocityPolicy(VelocityPolicy &&) = default;
  VelocityPolicy & operator=(VelocityPolicy &&) = delete;
  ~VelocityPolicy() = default;

  template<
    typename ...Args, int _rank = rank_,
    mpl::enable_if_t<_rank==1,int> = 0
    >
  VelocityPolicy(size_type extent0, Args && ...args)
    : projection_policy_t(std::forward<Args>(args)...),
      extents_{extent0}
  {}

  template<
    typename ...Args, int _rank = rank_,
    mpl::enable_if_t<_rank==2,int> = 0
    >
  VelocityPolicy(size_type extent0, size_type extent1, Args && ...args)
    : projection_policy_t(std::forward<Args>(args)...),
      extents_{extent0, extent1}
  {}

public:
  template <typename fom_system_t, int _rank = rank_>
  mpl::enable_if_t<_rank==1, galerkin_velocity_t>
  create(const fom_system_t & fomSystemObj) const
  {
    galerkin_velocity_t result(extents_[0]);
    ::pressio::ops::set_zero(result);
    return result;
  }

  template <typename fom_system_t, int _rank = rank_>
  mpl::enable_if_t<_rank==2, galerkin_velocity_t>
  create(const fom_system_t & fomSystemObj) const
  {
    galerkin_velocity_t result(extents_[0],extents_[1]);
    ::pressio::ops::set_zero(result);
    return result;
  }

  template<class galerkin_state_t, typename fom_system_t, typename scalar_t>
  void compute(const galerkin_state_t & galerkinState,
	       galerkin_velocity_t & galerkinRhs,
	       const fom_system_t  & fomSystemObj,
	       const scalar_t & timeToPassToRhsEvaluation) const
  {
    projection_policy_t::compute(galerkinRhs, galerkinState,
				 fomSystemObj, timeToPassToRhsEvaluation,
				 ::pressio::ode::n());
  }
};

}}}}//end namespace pressio::rom::galerkin::impl
#endif  // ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_POLICIES_ROM_GALERKIN_VELOCITY_POLICY_HPP_
