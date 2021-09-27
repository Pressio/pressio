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

#ifndef ROM_ROM_PUBLIC_API_LSPG_HPP_
#define ROM_ROM_PUBLIC_API_LSPG_HPP_

#include "./impl/rom_lspg_steady_compose.hpp"
#include "./impl/rom_lspg_unsteady_compose.hpp"

namespace pressio{ namespace rom{ namespace lspg{

// --------------------
// DEFAULT
// --------------------

// steady
template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType
  >
auto create_default_steady_problem(const FomSystemType & fomSysObj,
				   DecoderType & decoder,
				   const RomStateType & stateIn,
				   const FomReferenceStateType & fomRef)
{
  using return_t = typename impl::ComposeDefaultProblemSteady<
    FomSystemType, DecoderType, RomStateType, FomReferenceStateType
    >::type;

  return return_t(fomSysObj, decoder, stateIn, fomRef);
}

// steady, with preconditioner
template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class PreconditionerType
  >
auto create_default_steady_problem(const FomSystemType & fomSysObj,
				   DecoderType & decoder,
				   const RomStateType & stateIn,
				   const FomReferenceStateType & fomRef,
				   const PreconditionerType & prec)
{
  using return_t = typename impl::ComposePrecDefaultProblemSteady<
    FomSystemType, DecoderType, RomStateType, FomReferenceStateType, PreconditionerType
    >::type;

  return return_t(fomSysObj, decoder, stateIn, fomRef, prec);
}

// unsteady, cont-time
template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType
  >
auto create_default_unsteady_problem(::pressio::ode::StepScheme name,
				     const FomSystemType & fomSysObj,
				     DecoderType & decoder,
				     const RomStateType & stateIn,
				     const FomReferenceStateType & fomRef)
{
  using return_t = typename impl::ComposeDefaultProblemContTime<
    FomSystemType, DecoderType, RomStateType, FomReferenceStateType
    >::type;
  return return_t(name, fomSysObj, decoder, stateIn, fomRef);
}

// unsteady, cont-time, with preconditioner
template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class PreconditionerType
  >
auto create_default_unsteady_problem(::pressio::ode::StepScheme name,
				     const FomSystemType & fomSysObj,
				     DecoderType & decoder,
				     const RomStateType & stateIn,
				     const FomReferenceStateType & fomRef,
				     const PreconditionerType & prec)
{
  using return_t = typename impl::ComposePrecDefaultProblemContTime<
    FomSystemType, DecoderType, RomStateType,
    FomReferenceStateType, PreconditionerType
    >::type;

  return return_t(name, fomSysObj, decoder, stateIn, fomRef, prec);
}

// unsteady, discrete-time
template<
  std::size_t num_states,
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType
  >
auto create_default_unsteady_problem(const FomSystemType & fomSysObj,
				     DecoderType & decoder,
				     const RomStateType & stateIn,
				     const FomReferenceStateType & fomRef)
{
  using return_t = typename impl::ComposeDefaultProblemDiscTime<
    num_states, FomSystemType, DecoderType, RomStateType, FomReferenceStateType
    >::type;

  return return_t(::pressio::ode::StepScheme::ImplicitArbitrary,
		  fomSysObj, decoder, stateIn, fomRef);
}

// unsteady, discrete-time, with preconditioner
template<
  std::size_t num_states,
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class PreconditionerType
  >
auto create_default_unsteady_problem(const FomSystemType & fomSysObj,
				     DecoderType & decoder,
				     const RomStateType & stateIn,
				     const FomReferenceStateType & fomRef,
				     const PreconditionerType & prec)
{
  using return_t = typename impl::ComposePrecDefaultProblemDiscTime<
    num_states, FomSystemType, DecoderType, RomStateType,
    FomReferenceStateType, PreconditionerType>::type;

  return return_t(::pressio::ode::StepScheme::ImplicitArbitrary,
		  fomSysObj, decoder, stateIn, fomRef, prec);
}

// --------------------
// MASKED
// --------------------

// steady
template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class MaskerType
  >
auto create_masked_steady_problem(const FomSystemType & fomSysObj,
				  DecoderType & decoder,
				  const RomStateType & stateIn,
				  const FomReferenceStateType & fomRef,
				  const MaskerType & masker)
{

  using return_t = typename impl::ComposeMaskedProblemSteady<
    FomSystemType, DecoderType, RomStateType, FomReferenceStateType, MaskerType
    >::type;

  return return_t(fomSysObj, decoder, stateIn, fomRef, masker);
}

// steady, with preconditioner
template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class MaskerType,
  class PreconditionerType
  >
auto create_masked_steady_problem(const FomSystemType & fomSysObj,
				  DecoderType & decoder,
				  const RomStateType & stateIn,
				  const FomReferenceStateType & fomRef,
				  const MaskerType & masker,
				  const PreconditionerType & prec)
{

  using return_t = typename impl::ComposePrecMaskedProblemSteady<
    FomSystemType, DecoderType, RomStateType, FomReferenceStateType,
    MaskerType, PreconditionerType>::type;

  return return_t(fomSysObj, decoder, stateIn, fomRef, prec, masker);
}

// unsteady, cont-time
template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class MaskerType
  >
auto create_masked_unsteady_problem(::pressio::ode::StepScheme name,
				    const FomSystemType & fomSysObj,
				    DecoderType & decoder,
				    const RomStateType & stateIn,
				    const FomReferenceStateType & fomRef,
				    const MaskerType & masker)
{
  using return_t = typename impl::ComposeMaskedProblemContTime<
    FomSystemType, DecoderType, RomStateType, FomReferenceStateType, MaskerType
    >::type;

  return return_t(name, fomSysObj, decoder, stateIn, fomRef, masker);
}

// unsteady, cont-time, with preconditioner
template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class MaskerType,
  class PreconditionerType
  >
auto create_masked_unsteady_problem(::pressio::ode::StepScheme name,
				    const FomSystemType & fomSysObj,
				    DecoderType & decoder,
				    const RomStateType & stateIn,
				    const FomReferenceStateType & fomRef,
				    const MaskerType & masker,
				    const PreconditionerType & prec)
{
  using return_t = typename impl::ComposePrecMaskedProblemContTime<
    FomSystemType, DecoderType, RomStateType, FomReferenceStateType,
    MaskerType, PreconditionerType
    >::type;

  return return_t(name, fomSysObj, decoder, stateIn, fomRef, prec, masker);
}

// unsteady, discrete-time
template<
  std::size_t num_states,
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class MaskerType
  >
auto create_masked_unsteady_problem(const FomSystemType & fomSysObj,
				    DecoderType & decoder,
				    const RomStateType & stateIn,
				    const FomReferenceStateType & fomRef,
				    const MaskerType & masker)
{

  using return_t = typename impl::ComposeMaskedProblemDiscTime<
    num_states, FomSystemType, DecoderType, RomStateType,
    FomReferenceStateType, MaskerType>::type;

  return return_t(::pressio::ode::StepScheme::ImplicitArbitrary,
		  fomSysObj, decoder, stateIn, fomRef, masker);
}

// unsteady, discrete-time, with preconditioner
template<
  std::size_t num_states,
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class MaskerType,
  class PreconditionerType
  >
auto create_masked_unsteady_problem(const FomSystemType & fomSysObj,
				    DecoderType & decoder,
				    const RomStateType & stateIn,
				    const FomReferenceStateType & fomRef,
				    const MaskerType & masker,
				    const PreconditionerType & prec)
{

  using return_t = typename impl::ComposePrecMaskedProblemDiscTime<
    num_states, FomSystemType, DecoderType, RomStateType, FomReferenceStateType,
    MaskerType, PreconditionerType>::type;

  return return_t(::pressio::ode::StepScheme::ImplicitArbitrary,
		  fomSysObj, decoder, stateIn, fomRef, prec, masker);
}

// --------------------
// HYPRED
// --------------------

// steady
template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType
  >
auto create_hyperreduced_steady_problem(const FomSystemType & fomSysObj,
					DecoderType & decoder,
					const RomStateType & stateIn,
					const FomReferenceStateType & fomRef)
{

  using return_t = typename impl::ComposeHyperreducedProblemSteady<
    FomSystemType, DecoderType, RomStateType, FomReferenceStateType
    >::type;

  return return_t(fomSysObj, decoder, stateIn, fomRef);
}

// steady, with precond
template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class PreconditionerType
  >
auto create_hyperreduced_steady_problem(const FomSystemType & fomSysObj,
					DecoderType & decoder,
					const RomStateType & stateIn,
					const FomReferenceStateType & fomRef,
					const PreconditionerType & prec)
{
  using ReturnType = typename impl::ComposePrecHypredProblemSteady<
    FomSystemType, DecoderType, RomStateType, FomReferenceStateType, PreconditionerType
    >::type;

  return ReturnType(fomSysObj, decoder, stateIn, fomRef, prec);
}

// unsteady, cont-time
template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class HypRedOperatorUpdaterType
  >
auto create_hyperreduced_unsteady_problem(::pressio::ode::StepScheme name,
					  const FomSystemType & fomSysObj,
					  DecoderType & decoder,
					  const RomStateType & stateIn,
					  const FomReferenceStateType & fomRef,
					  const HypRedOperatorUpdaterType & operatorUpdater)
{

  using return_t = typename impl::ComposeHypRedProblemContTime<
    FomSystemType, DecoderType, RomStateType,
    FomReferenceStateType, HypRedOperatorUpdaterType>::type;

  return return_t(name, fomSysObj, decoder, stateIn, fomRef, operatorUpdater);
}

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
// unsteady, cont-time
// this overload is specific to trilinos for now, where we can use the underlying
// maps of the operators to figure out how to combine them.
// if you call this for non-trilinos data types, you get an error
template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType
  >
auto create_hyperreduced_unsteady_problem(::pressio::ode::StepScheme name,
					  const FomSystemType & fomSysObj,
					  DecoderType & decoder,
					  const RomStateType & stateIn,
					  const FomReferenceStateType & fomRef)
{

  // this is the cont-time, we can query the adapter's types and
  // error out if they are not trilinos types
  using fom_state_type = typename mpl::remove_cvref_t<FomSystemType>::state_type;
  using fom_velo_type = typename mpl::remove_cvref_t<FomSystemType>::velocity_type;
  static_assert
    (pressio::Traits<fom_state_type>::package_identifier ==
     pressio::PackageIdentifier::Trilinos and
     pressio::Traits<fom_velo_type>::package_identifier ==
     pressio::PackageIdentifier::Trilinos,
     "You are calling an overload of create_hyperreduced_unsteaady_problem that \
is only valid for Trilinos type");

  using hr_operator_updater_t = impl::HypRedUpdaterTrilinos;
  using return_t = typename impl::ComposeHypRedProblemContTime<
    FomSystemType, DecoderType, RomStateType,
    FomReferenceStateType, hr_operator_updater_t>::type;

  return return_t(name, fomSysObj, decoder, stateIn,
		  fomRef, hr_operator_updater_t());
}

// unsteady, cont-time
// this overload is specific to trilinos for now, where we can use the underlying
// maps of the operators to figure out how to combine them.
// if you call this for non-trilinos data types, you get an error
template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class PreconditionerType
  >
auto create_prec_hyperreduced_unsteady_problem(::pressio::ode::StepScheme name,
					       const FomSystemType & fomSysObj,
					       DecoderType & decoder,
					       const RomStateType & stateIn,
					       const FomReferenceStateType & fomRef,
					       const PreconditionerType & prec)
{

  // this is the cont-time, we can query the adapter's types and
  // error out if they are not trilinos types
  using fom_state_type = typename mpl::remove_cvref_t<FomSystemType>::state_type;
  using fom_velo_type = typename mpl::remove_cvref_t<FomSystemType>::velocity_type;
  static_assert
    (pressio::Traits<fom_state_type>::package_identifier ==
     pressio::PackageIdentifier::Trilinos and
     pressio::Traits<fom_velo_type>::package_identifier ==
     pressio::PackageIdentifier::Trilinos,
     "You are calling an overload of create_hyperreduced_unsteaady_problem that \
is only valid for Trilinos type");

  using hr_operator_updater_t = impl::HypRedUpdaterTrilinos;
  using return_t = typename impl::ComposePrecHypRedProblemContTime<
    FomSystemType, DecoderType, RomStateType, FomReferenceStateType,
    hr_operator_updater_t, PreconditionerType>::type;

  return return_t(name, fomSysObj, decoder, stateIn, fomRef,
		  hr_operator_updater_t(), prec);
}
#endif

// unsteady, cont-time, with preconditioner
template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class HypRedOperatorUpdaterType,
  class PreconditionerType
  >
auto create_prec_hyperreduced_unsteady_problem(::pressio::ode::StepScheme name,
					       const FomSystemType & fomSysObj,
					       DecoderType & decoder,
					       const RomStateType & stateIn,
					       const FomReferenceStateType & fomRef,
					       const HypRedOperatorUpdaterType & operatorUpdater,
					       const PreconditionerType & prec)
{

  using return_t = typename impl::ComposePrecHypRedProblemContTime<
    FomSystemType, DecoderType, RomStateType, FomReferenceStateType,
    HypRedOperatorUpdaterType, PreconditionerType>::type;

  return return_t(name, fomSysObj, decoder, stateIn, fomRef, operatorUpdater, prec);
}

// unsteady, discrete-time
template<
  std::size_t num_states,
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType
  >
auto create_hyperreduced_unsteady_problem(const FomSystemType & fomSysObj,
					  DecoderType & decoder,
					  const RomStateType & stateIn,
					  const FomReferenceStateType & fomRef)
{

  using return_t = typename impl::ComposeHypRedProblemDiscTime<
    num_states, FomSystemType, DecoderType, RomStateType, FomReferenceStateType
    >::type;

  return return_t(::pressio::ode::StepScheme::ImplicitArbitrary,
		  fomSysObj, decoder, stateIn, fomRef);
}

// unsteady, discrete-time, with preconditioner
template<
  std::size_t num_states,
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class PreconditionerType
  >
auto create_prec_hyperreduced_unsteady_problem(const FomSystemType & fomSysObj,
					  DecoderType & decoder,
					  const RomStateType & stateIn,
					  const FomReferenceStateType & fomRef,
					  const PreconditionerType & prec)
{

  using return_t = typename impl::ComposePrecHypRedProblemDiscTime<
    num_states, FomSystemType, DecoderType, RomStateType,
    FomReferenceStateType, PreconditionerType>::type;

  return return_t(::pressio::ode::StepScheme::ImplicitArbitrary,
		  fomSysObj, decoder, stateIn, fomRef, prec);
}

}}}//end namespace pressio::rom::lspg
#endif  // ROM_ROM_PUBLIC_API_LSPG_HPP_
