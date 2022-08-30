/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_cont_time_decorators.hpp
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

#ifndef ROM_IMPL_ROM_LSPG_UNSTEADY_CONT_TIME_DECORATORS_HPP_
#define ROM_IMPL_ROM_LSPG_UNSTEADY_CONT_TIME_DECORATORS_HPP_

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{

template <class PreconditionerType, class PreconditionableType>
class PrecDecoratorResidual : public PreconditionableType
{
public:
  // required
  using residual_type = typename PreconditionableType::residual_type;

public:
  PrecDecoratorResidual() = delete;
  PrecDecoratorResidual(const PrecDecoratorResidual &) = default;
  PrecDecoratorResidual & operator=(const PrecDecoratorResidual &) = default;
  PrecDecoratorResidual(PrecDecoratorResidual &&) = default;
  PrecDecoratorResidual & operator=(PrecDecoratorResidual &&) = default;
  ~PrecDecoratorResidual() = default;

  template <class ... Args>
  PrecDecoratorResidual(
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
			pybind11::object preconditioner,
#else
			const PreconditionerType & preconditioner,
#endif
			Args && ... args)
    : PreconditionableType(std::forward<Args>(args)...),
      preconditionerObj_(preconditioner)
  {}

public:
  residual_type create() const{ return PreconditionableType::create(); }

  template <
    class LspgStateType,
    class LspgStencilStatesContainerType,
    class LspgStencilVelocitiesContainerType,
    class ScalarType,
    class StepType
    >
  void operator()(::pressio::ode::StepScheme name,
	       const LspgStateType & lspgState,
	       const LspgStencilStatesContainerType & lspgStencilStates,
	       LspgStencilVelocitiesContainerType & lspgStencilVelocities,
	       const ScalarType & t_np1,
	       const ScalarType & dt,
	       const StepType & currentStepNumber,
	       residual_type & unprecLspgResidual) const
  {
    PreconditionableType::operator()(name, lspgState, lspgStencilStates,
				     lspgStencilVelocities,
				     t_np1, dt, currentStepNumber, unprecLspgResidual);
    const auto & fomState = fomStatesMngr_(::pressio::ode::nPlusOne());
    preconditionerObj_(fomState, t_np1, unprecLspgResidual);
  }

private:
  using PreconditionableType::fomStatesMngr_;
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  const PreconditionerType preconditionerObj_;
#else
  std::reference_wrapper<const PreconditionerType> preconditionerObj_;
#endif
};


template <class PreconditionerType, class PreconditionableType>
class PrecDecoratorJacobian : public PreconditionableType
{
public:
  // required
  using jacobian_type = typename PreconditionableType::jacobian_type;

public:
  PrecDecoratorJacobian() = delete;
  PrecDecoratorJacobian(const PrecDecoratorJacobian &) = default;
  PrecDecoratorJacobian & operator=(const PrecDecoratorJacobian &) = default;
  PrecDecoratorJacobian(PrecDecoratorJacobian &&) = default;
  PrecDecoratorJacobian & operator=(PrecDecoratorJacobian &&) = default;
  ~PrecDecoratorJacobian() = default;

  template <class ... Args>
  PrecDecoratorJacobian(
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
			pybind11::object preconditioner,
#else
			const PreconditionerType & preconditioner,
#endif
			Args && ... args)
    : PreconditionableType(std::forward<Args>(args)...),
      preconditionerObj_(preconditioner)
  {}

public:
  jacobian_type create() const{ return PreconditionableType::create(); }

  template <
    class LspgStateType,
    class LspgStencilStatesContainerType,
    class ScalarType,
    class StepType
    >
  void operator()(::pressio::ode::StepScheme name,
	       const LspgStateType & lspgState,
	       const LspgStencilStatesContainerType & lspgStencilStates,
	       const ScalarType & t_np1,
	       const ScalarType & dt,
	       const StepType & currentStepNumber,
	       jacobian_type & unprecLspgJacobian) const
  {
    PreconditionableType::operator()(name, lspgState, lspgStencilStates,
				     t_np1, dt, currentStepNumber, unprecLspgJacobian);
    const auto & fomState = fomStatesMngr_(::pressio::ode::nPlusOne());
    preconditionerObj_(fomState, t_np1, unprecLspgJacobian);
  }

private:
  using PreconditionableType::fomStatesMngr_;
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  const PreconditionerType preconditionerObj_;
#else
  std::reference_wrapper<const PreconditionerType> preconditionerObj_;
#endif
};


template <class MaskerType, class MaskableType>
class MaskDecoratorResidual : public MaskableType
{
public:
  // required
  using residual_type = typename MaskableType::residual_type;

public:
  MaskDecoratorResidual() = delete;
  MaskDecoratorResidual(const MaskDecoratorResidual &) = default;
  MaskDecoratorResidual & operator=(const MaskDecoratorResidual &) = default;
  MaskDecoratorResidual(MaskDecoratorResidual &&) = default;
  MaskDecoratorResidual & operator=(MaskDecoratorResidual &&) = default;
  ~MaskDecoratorResidual() = default;

  template <class ... Args>
  MaskDecoratorResidual(
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
			pybind11::object masker,
#else
			const MaskerType & masker,
#endif
			Args && ... args)
    : MaskableType(std::forward<Args>(args)...),
      unmaskedResidual_(MaskableType::create()),
      masker_(masker)
  {}

public:
  residual_type create() const{
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
    residual_type R(masker_.createApplyMaskResult(unmaskedResidual_));
#else
    residual_type R(masker_.get().createApplyMaskResult(unmaskedResidual_));
#endif
    return R;
  }

  template <
    class LspgStateType,
    class LspgStencilStatesContainerType,
    class LspgStencilVelocitiesContainerType,
    class ScalarType,
    class StepType
    >
  void operator()(::pressio::ode::StepScheme name,
	       const LspgStateType & lspgState,
	       const LspgStencilStatesContainerType & lspgStencilStates,
	       LspgStencilVelocitiesContainerType & lspgStencilVelocities,
	       const ScalarType & t_np1,
	       const ScalarType & dt,
	       const StepType & currentStepNumber,
	       residual_type & lspgResidual) const
  {
    MaskableType::operator()(name, lspgState, lspgStencilStates, lspgStencilVelocities,
			  t_np1, dt, currentStepNumber, unmaskedResidual_);
    masker_(unmaskedResidual_, t_np1, lspgResidual);
  }

private:
  mutable residual_type unmaskedResidual_;

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  const MaskerType masker_;
#else
  std::reference_wrapper<const MaskerType> masker_;
#endif
};

template <class MaskerType, class MaskableType>
class MaskDecoratorJacobian : public MaskableType
{
public:
  // required
  using jacobian_type = typename MaskableType::jacobian_type;

public:
  MaskDecoratorJacobian() = delete;
  MaskDecoratorJacobian(const MaskDecoratorJacobian &) = default;
  MaskDecoratorJacobian & operator=(const MaskDecoratorJacobian &) = default;
  MaskDecoratorJacobian(MaskDecoratorJacobian &&) = default;
  MaskDecoratorJacobian & operator=(MaskDecoratorJacobian &&) = default;
  ~MaskDecoratorJacobian() = default;

  template <class ... Args>
  MaskDecoratorJacobian(
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
			pybind11::object masker,
#else
			const MaskerType & masker,
#endif
			Args && ... args)
    : MaskableType(std::forward<Args>(args)...),
      unmaskedJacobian_(MaskableType::create()),
      masker_(masker)
  {}

public:
  jacobian_type create() const{
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
    jacobian_type R(masker_.createApplyMaskResult(unmaskedJacobian_));
#else
    jacobian_type R(masker_.get().createApplyMaskResult(unmaskedJacobian_));
#endif
    return R;
  }

  template <
    class LspgStateType,
    class LspgStencilStatesContainerType,
    class ScalarType,
    class StepType
    >
  void operator()(::pressio::ode::StepScheme name,
	       const LspgStateType & lspgState,
	       const LspgStencilStatesContainerType & lspgStencilStates,
	       const ScalarType & t_np1,
	       const ScalarType & dt,
	       const StepType & currentStepNumber,
	       jacobian_type & lspgJacobian) const
  {
    MaskableType::operator()(name, lspgState, lspgStencilStates,
			  t_np1, dt, currentStepNumber, unmaskedJacobian_);
    masker_(unmaskedJacobian_, t_np1, lspgJacobian);
  }

private:
  mutable jacobian_type unmaskedJacobian_;
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  const MaskerType masker_;
#else
  std::reference_wrapper<const MaskerType> masker_;
#endif
};

}}}} //end namespace pressio::rom::lspg::decorator
#endif  // ROM_IMPL_ROM_LSPG_UNSTEADY_CONT_TIME_DECORATORS_HPP_
