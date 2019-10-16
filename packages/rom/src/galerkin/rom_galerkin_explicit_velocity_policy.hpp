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

#ifndef ROM_DEFAULT_GALERKIN_EXPLICIT_VELOCITY_POLICY_HPP_
#define ROM_DEFAULT_GALERKIN_EXPLICIT_VELOCITY_POLICY_HPP_

#include "../rom_fwd.hpp"
#include "../../../ode/src/explicit/policies/ode_explicit_velocity_policy_base.hpp"
#include "../rom_data_fom_states.hpp"

namespace pressio{ namespace rom{

template <
  typename fom_states_data_type,
  typename fom_rhs_t,
  typename decoder_t,
  typename ud_ops
  >
class DefaultGalerkinExplicitVelocityPolicy
  : public ode::policy::ExplicitVelocityPolicyBase<
       DefaultGalerkinExplicitVelocityPolicy<fom_states_data_type,
					     fom_rhs_t,
					     decoder_t, ud_ops>>
{

protected:
  using this_t = DefaultGalerkinExplicitVelocityPolicy<fom_states_data_type,
						       fom_rhs_t,
						       decoder_t,
						       ud_ops>;
  friend ode::policy::ExplicitVelocityPolicyBase<this_t>;

public:
  static constexpr bool isResidualPolicy_ = true;

public:
  DefaultGalerkinExplicitVelocityPolicy() = delete;
  ~DefaultGalerkinExplicitVelocityPolicy() = default;

  DefaultGalerkinExplicitVelocityPolicy(const fom_rhs_t & fomRhs,
					fom_states_data_type & fomStates,
					const decoder_t & decoder)
    : R_{fomRhs},
      fomStates_(fomStates),
      decoder_(decoder){}

// #ifdef PRESSIO_ENABLE_TPL_PYBIND11
//   // this cnstr only enabled when udOps is non-void and python
//   template <
//     typename _ud_ops = ud_ops,
//     mpl::enable_if_t<
//       mpl::is_same<_ud_ops, pybind11::object>::value
//       > * = nullptr
//     >
//   DefaultGalerkinExplicitVelocityPolicy(fom_states_data_type & fomStates,
// 					const fom_rhs_t & fomResids,
// 					const decoder_t & decoder,
// 					const _ud_ops & udOps)
//     : fomStates_(fomStates),
//       fom_rhs_t(fomResids),
//       decoder_(decoder),
//       udOps_{udOps}{}
// #endif

public:
  /* for now, the ROM state and ROM velocity must be of the same type */
  template <
    typename galerkin_state_t,
    typename fom_t,
    typename scalar_t
  >
  void operator()(const galerkin_state_t  & romY,
		  galerkin_state_t	  & romR,
  		  const fom_t		  & app,
		  scalar_t		  t) const{
    this->compute_impl(romY, romR, app, t);
  }

  template <
    typename galerkin_state_t,
    typename fom_t,
    typename scalar_t
    >
  galerkin_state_t operator()(const galerkin_state_t  & romY,
			      const fom_t	      & app,
			      scalar_t		 t) const
  {
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
    // TODO: make this better, maybe initialized somewhere else
    galerkin_state_t result( const_cast<galerkin_state_t &>(romY).request() );
#else
    galerkin_state_t result(romY);
#endif
    ::pressio::containers::ops::set_zero(result);
    this->compute_impl(romY, result, app, t);
    return result;
  }


private:

  template <
  typename galerkin_state_t,
  typename fom_t,
  typename scalar_t,
  typename _ud_ops = ud_ops,
  mpl::enable_if_t<
    std::is_void<_ud_ops>::value and
    ::pressio::containers::meta::is_vector_wrapper<galerkin_state_t>::value
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
    and !::pressio::containers::meta::is_array_pybind11<galerkin_state_t>::value
#endif
    > * = nullptr
  >
  void compute_impl(const galerkin_state_t  & romY,
		    galerkin_state_t	    & romR,
		    const fom_t		    & app,
		    scalar_t		    t) const
  {
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("galerkin explicit velocity");
#endif

    fomStates_.template reconstructCurrentFomState(romY);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->start("fom eval rhs");
#endif
    const auto & yFom = fomStates_.getCRefToFomState();
    app.velocity(*yFom.data(), t, *R_.data());

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("fom eval rhs");
    timer->start("phiT*fomRhs");
#endif

    const auto & phi = decoder_.getReferenceToJacobian();
    containers::ops::dot(phi, R_, romR);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("phiT*fomRhs");
    timer->stop("galerkin explicit velocity");
#endif
  }


#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  template <
  typename galerkin_state_t,
  typename fom_t,
  typename scalar_t,
  typename _ud_ops = ud_ops,
  mpl::enable_if_t<
    ::pressio::mpl::is_same<_ud_ops, pybind11::object>::value and
    ::pressio::containers::meta::is_array_pybind11<galerkin_state_t>::value
    > * = nullptr
  >
  void compute_impl(const galerkin_state_t  & romY,
		    galerkin_state_t	    & romR,
		    const fom_t		    & app,
		    scalar_t		    t) const
  {
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("galerkin explicit velocity");
#endif

    fomStates_.template reconstructCurrentFomState(romY);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->start("fom eval rhs");
#endif
    app.attr("velocity")(yFom_, t, fomRhs_);
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("fom eval rhs");
    timer->start("phiT*fomRhs");
#endif
    const auto & phi = decoder_.getReferenceToJacobian();
    auto constexpr transposePhi = true;
    udOps_.attr("multiply2")(phi, fomRhs_, romR, transposePhi);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("phiT*fomRhs");
    timer->stop("galerkin explicit velocity");
#endif
  }
#endif


protected:
  mutable fom_rhs_t R_ = {};
  const decoder_t & decoder_;
  fom_states_data_type & fomStates_;

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  typename std::conditional<
    mpl::is_same<ud_ops, pybind11::object>::value,
    ud_ops,
    const ud_ops *
    >::type udOps_ = {};
#else
    const ud_ops * udOps_ = {};
#endif


};//end class

}}//end namespace pressio::rom
#endif
