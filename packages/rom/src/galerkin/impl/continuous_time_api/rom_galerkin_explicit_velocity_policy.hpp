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

#ifndef ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_ROM_GALERKIN_EXPLICIT_VELOCITY_POLICY_HPP_
#define ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_ROM_GALERKIN_EXPLICIT_VELOCITY_POLICY_HPP_

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

template <
  typename galerkin_state_t,
  typename fom_states_manager_t,
  typename fom_rhs_t,
  typename decoder_t,
  typename ud_ops
  >
class ExplicitVelocityPolicy
{

public:
  static constexpr bool isResidualPolicy_ = true;

public:
  ExplicitVelocityPolicy() = delete;
  ExplicitVelocityPolicy(const ExplicitVelocityPolicy &) = default;
  ExplicitVelocityPolicy & operator=(const ExplicitVelocityPolicy &) = default;
  ExplicitVelocityPolicy(ExplicitVelocityPolicy &&) = default;
  ExplicitVelocityPolicy & operator=(ExplicitVelocityPolicy &&) = default;
  ~ExplicitVelocityPolicy() = default;

  // 1. void ops
  template <
    typename _fom_rhs_t = fom_rhs_t,
    typename _ud_ops = ud_ops,
    ::pressio::mpl::enable_if_t<
      std::is_void<_ud_ops>::value, int
      > = 0
    >
  ExplicitVelocityPolicy(const _fom_rhs_t & fomRhs,
			 fom_states_manager_t & fomStatesMngr,
			 const decoder_t & decoder)
    : fomRhs_{fomRhs},
      decoder_{decoder},
      phi_(decoder.jacobianCRef()),
      fomStatesMngr_(fomStatesMngr){}

  // 2. non-void ops
  template <
    typename _fom_rhs_t = fom_rhs_t,
    typename _ud_ops = ud_ops,
    ::pressio::mpl::enable_if_t<
      !std::is_void<_ud_ops>::value, int
      > = 0
    >
  ExplicitVelocityPolicy(const _fom_rhs_t & fomRhs,
			 fom_states_manager_t & fomStatesMngr,
			 const decoder_t & decoder,
			 const _ud_ops & udOps)
    : fomRhs_{fomRhs},
      decoder_{decoder},
      phi_(decoder.jacobianCRef()),
      fomStatesMngr_(fomStatesMngr),
      udOps_{&udOps}{}

public:
  template <typename fom_t>
  galerkin_state_t create(const fom_t & app) const
  {
    galerkin_state_t result(phi_.get().extent(1));
    ::pressio::ops::set_zero(result);
    return result;
  }

  /* for now, the ROM state and ROM velocity must be of the same type */
  template <typename fom_t, typename scalar_t>
  void compute(const galerkin_state_t & romState,
	       galerkin_state_t & romRhs,
	       const fom_t	& app,
	       const scalar_t & t) const
  {
    this->compute_impl(romState, romRhs, app, t);
  }

private:
  // --------------------------------------------
  // compute RHS of ode
  // --------------------------------------------
  template <
    typename scalar_t,
    typename result_t,
    typename _ops_t = ud_ops
  >
  ::pressio::mpl::enable_if_t< std::is_void<_ops_t>::value >
  applyDecoderJacobianToFomVel(result_t & result) const
  {
    constexpr auto zero = ::pressio::utils::constants<scalar_t>::zero();
    constexpr auto one  = ::pressio::utils::constants<scalar_t>::one();
    ::pressio::ops::product(::pressio::transpose(), one, phi_.get(),
			    fomRhs_, zero, result);
  }

  template <
    typename scalar_t,
    typename result_t,
    typename _ops_t = ud_ops
  >
  ::pressio::mpl::enable_if_t< !std::is_void<_ops_t>::value >
  applyDecoderJacobianToFomVel(result_t & result) const
  {
    constexpr auto zero = ::pressio::utils::constants<scalar_t>::zero();
    constexpr auto one  = ::pressio::utils::constants<scalar_t>::one();
    udOps_->product(::pressio::transpose(), one, *(phi_.get().data()),
		    *fomRhs_.data(), zero, result);
  }

  template <typename fom_system_t, typename scalar_t>
  void compute_impl(const galerkin_state_t  & romState,
		    galerkin_state_t	    & romRhs,
		    const fom_system_t	    & fomSystemObj,
		    const scalar_t	    & time) const
  {
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("galerkin explicit velocity");
#endif

    // any time compute_impl is called, it means the romState
    // has changed, so tell decoder to update the Jacobian
    decoder_.get().updateJacobian(romState);

    // reconstruct the current fom state
    fomStatesMngr_.get().reconstructCurrentFomState(romState);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->start("fom eval rhs");
#endif
    const auto & yFom = fomStatesMngr_.get().currentFomStateCRef();
    fomSystemObj.velocity(*yFom.data(), time, *fomRhs_.data());

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("fom eval rhs");
    timer->start("phiT*fomRhs");
#endif

    // apply decoder's jacobian to velocity
    (*this).template applyDecoderJacobianToFomVel<scalar_t>(romRhs);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("phiT*fomRhs");
    timer->stop("galerkin explicit velocity");
#endif
  }

protected:
  mutable fom_rhs_t fomRhs_ = {};
  std::reference_wrapper<const decoder_t> decoder_;
  std::reference_wrapper<const typename decoder_t::jacobian_type> phi_;
  std::reference_wrapper<fom_states_manager_t> fomStatesMngr_;
  const ud_ops * udOps_ = {};

};//end class

}}}}//end namespace pressio::rom::galerkin::impl
#endif  // ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_ROM_GALERKIN_EXPLICIT_VELOCITY_POLICY_HPP_
