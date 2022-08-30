/*
//@HEADER
// ************************************************************************
//
// rom_public_api_lspg.hpp
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

#ifndef ROM_ROM_IMPL_RESIDUAL_EVALUATOR_HPP_
#define ROM_ROM_IMPL_RESIDUAL_EVALUATOR_HPP_

namespace pressio{ namespace rom{ namespace impl{

template <class fom_state_reconstr_type>
class ResidualEvaluator
{

  using fom_state_type = typename fom_state_reconstr_type::fom_state_type;
  ::pressio::ode::StepScheme step_scheme_;
  const fom_state_reconstr_type fomRec_;

public:
  ResidualEvaluator() = delete;

  ResidualEvaluator(::pressio::ode::StepScheme step_scheme,
		    const typename fom_state_reconstr_type::decoder_type & decoder,
		    const fom_state_type & fomRef)
    : step_scheme_(step_scheme), fomRec_(fomRef, decoder)
  {}

  template<class FomAppType, class RomStatesProducer, class ResidualType, class ObserverType>
  void operator()(const FomAppType & fomObj,
		  const RomStatesProducer & romStatesProducer,
		  ResidualType & fomResidual,
		  ObserverType && observer) const
  {
    if (step_scheme_ == ::pressio::ode::StepScheme::BDF1){
      return this->bdf1impl(fomObj, romStatesProducer, fomResidual,
			    std::forward<ObserverType>(observer));
    }

    else if (step_scheme_ == ::pressio::ode::StepScheme::BDF2){
      return this->bdf2impl(fomObj, romStatesProducer, fomResidual,
			    std::forward<ObserverType>(observer));
    }

    else if (step_scheme_ == ::pressio::ode::StepScheme::CrankNicolson){
      return this->cnimpl(fomObj, romStatesProducer, fomResidual,
			  std::forward<ObserverType>(observer));
    }

    else{
      throw std::runtime_error("Invalid scheme");
    }
  }

private:
  template<class FomAppType, class RomStatesProducer, class ResidualType, class ObserverType>
  void bdf1impl(const FomAppType & fomObj,
		const RomStatesProducer & romStatesProducer,
		ResidualType & fomResidual,
		ObserverType && observer) const
  {
    ::pressio::ops::set_zero(fomResidual);
    const std::size_t numInstances = romStatesProducer.numberOfTimeSchemeStencilInstances();

    for (std::size_t instanceIt=0; instanceIt<numInstances; ++instanceIt){
      //
      // r(y_n) = y_n - y_n-1 - dt F(y_n, t_n)
      //
      auto romState_nm1 = romStatesProducer.stateAt(instanceIt, ::pressio::ode::nMinusOne());
      auto fomState_nm1 = fomRec_(romState_nm1);

      auto romState_n = romStatesProducer.stateAt(instanceIt, ::pressio::ode::n());
      auto fomState_n = fomRec_(romState_n);

      const auto time_n = romStatesProducer.timeAt(instanceIt, ::pressio::ode::n());
      const auto dt   = romStatesProducer.timeStepSizeAt(instanceIt);

      fomObj.velocity(fomState_n, time_n, fomResidual);

      using scalar_type = typename ::pressio::Traits<ResidualType>::scalar_type;
      constexpr auto one = static_cast<scalar_type>(1);
      ::pressio::ops::update(fomResidual, -dt,
			     fomState_n, one,
			     fomState_nm1, -one);

      observer(instanceIt, fomResidual);
    }
  }

  template<class FomAppType, class RomStatesProducer, class ResidualType, class ObserverType>
  void bdf2impl(const FomAppType & fomObj,
		const RomStatesProducer & romStatesProducer,
		ResidualType & fomResidual,
		ObserverType && observer) const
  {
    ::pressio::ops::set_zero(fomResidual);
    const std::size_t numInstances = romStatesProducer.numberOfTimeSchemeStencilInstances();

    for (std::size_t instanceIt=0; instanceIt<numInstances; ++instanceIt){
      //
      // r(y_n) = y_n - (4/3) y_n-1 + (1/3) y_n-2 - dt F(y_n, t_n)
      //
      auto romState_nm2 = romStatesProducer.stateAt(instanceIt, ::pressio::ode::nMinusTwo());
      auto fomState_nm2 = fomRec_(romState_nm2);

      auto romState_nm1 = romStatesProducer.stateAt(instanceIt, ::pressio::ode::nMinusOne());
      auto fomState_nm1 = fomRec_(romState_nm1);

      auto romState_n = romStatesProducer.stateAt(instanceIt, ::pressio::ode::n());
      auto fomState_n = fomRec_(romState_n);

      const auto time_n = romStatesProducer.timeAt(instanceIt, ::pressio::ode::n());
      const auto dt   = romStatesProducer.timeStepSizeAt(instanceIt);
      fomObj.velocity(fomState_n, time_n, fomResidual);

      using scalar_type = typename ::pressio::Traits<ResidualType>::scalar_type;
      constexpr auto one = static_cast<scalar_type>(1);
      constexpr auto two = static_cast<scalar_type>(2);
      constexpr auto three = static_cast<scalar_type>(3);
      constexpr auto four  = two*two;

      ::pressio::ops::update(fomResidual, -(two/three)*dt,
			     fomState_n, one,
			     fomState_nm1, -four/three,
			     fomState_nm2,  one/three);

      observer(instanceIt, fomResidual);
    }
  }

  template<class FomAppType, class RomStatesProducer, class ResidualType, class ObserverType>
  void cnimpl(const FomAppType & fomObj,
	      const RomStatesProducer & romStatesProducer,
	      ResidualType & fomResidual,
	      ObserverType && observer) const
  {
    auto fomVel = ::pressio::ops::clone(fomResidual);
    ::pressio::ops::set_zero(fomResidual);
    ::pressio::ops::set_zero(fomVel);

    const std::size_t numInstances = romStatesProducer.numberOfTimeSchemeStencilInstances();
    for (std::size_t instanceIt=0; instanceIt<numInstances; ++instanceIt){
      //
      // r(y_n) = y_n - y_n-1 - 0.5 dt (F(y_n, t_n) + F(y_nm1, t_nm1))
      //
      auto romState_nm1 = romStatesProducer.stateAt(instanceIt, ::pressio::ode::nMinusOne());
      auto fomState_nm1 = fomRec_(romState_nm1);

      auto romState_n = romStatesProducer.stateAt(instanceIt, ::pressio::ode::n());
      auto fomState_n = fomRec_(romState_n);

      const auto time_nm1 = romStatesProducer.timeAt(instanceIt, ::pressio::ode::nMinusOne());
      const auto time_n   = romStatesProducer.timeAt(instanceIt, ::pressio::ode::n());
      const auto dt   = romStatesProducer.timeStepSizeAt(instanceIt);

      fomObj.velocity(fomState_n, time_n, fomResidual);
      fomObj.velocity(fomState_nm1, time_nm1, fomVel);

      using scalar_type = typename ::pressio::Traits<ResidualType>::scalar_type;
      constexpr auto one = static_cast<scalar_type>(1);
      constexpr auto two = static_cast<scalar_type>(2);
      ::pressio::ops::update(fomResidual, -(one/two)*dt,
			     fomVel, -(one/two)*dt,
			     fomState_n, one,
			     fomState_nm1, -one);

      observer(instanceIt, fomResidual);
    }
  }

};

}}} // end namespace pressio::rom::impl
#endif
