/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_discrete_time_decorators.hpp
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

#ifndef ROM_IMPL_ROM_LSPG_UNSTEADY_DISCRETE_TIME_DECORATORS_HPP_
#define ROM_IMPL_ROM_LSPG_UNSTEADY_DISCRETE_TIME_DECORATORS_HPP_

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{

template <class PreconditionerType, class preconditionable>
struct PrecDecoratorDiscreteTimeSystem : preconditionable
{
  using scalar_type = typename preconditionable::scalar_type;
  using state_type  = typename preconditionable::state_type;
  using discrete_time_residual_type = typename preconditionable::discrete_time_residual_type;
  using discrete_time_jacobian_type = typename preconditionable::discrete_time_jacobian_type;

  PrecDecoratorDiscreteTimeSystem() = delete;
  PrecDecoratorDiscreteTimeSystem(const PrecDecoratorDiscreteTimeSystem &) = default;
  PrecDecoratorDiscreteTimeSystem & operator=(const PrecDecoratorDiscreteTimeSystem &) = delete;
  PrecDecoratorDiscreteTimeSystem(PrecDecoratorDiscreteTimeSystem &&) = default;
  PrecDecoratorDiscreteTimeSystem & operator=(PrecDecoratorDiscreteTimeSystem &&) = delete;
  ~PrecDecoratorDiscreteTimeSystem() = default;

  template <class ... Args>
  PrecDecoratorDiscreteTimeSystem(
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
			pybind11::object preconditioner,
#else
			const PreconditionerType & preconditioner,
#endif
			Args && ... args)
    : preconditionable(std::forward<Args>(args)...),
      precFunctor_(preconditioner)
  {}

  discrete_time_residual_type createDiscreteTimeResidual() const{
    return preconditionable::createDiscreteTimeResidual();
  }

  discrete_time_jacobian_type createDiscreteTimeJacobian() const{
    return preconditionable::createDiscreteTimeJacobian();
  }

  template<class StepCountType, class ...Args>
  void discreteTimeResidual(StepCountType currentStepNumber,
                            scalar_type time_np1,
                            scalar_type dt,
                            discrete_time_residual_type & R,
			    Args && ...args) const
  {
    preconditionable::discreteTimeResidual(currentStepNumber, time_np1, dt,
					   R, std::forward<Args>(args)...);
    const auto & ynp1 = fomStatesMngr_(::pressio::ode::nPlusOne());
    precFunctor_(ynp1, time_np1, R);
  }

  template<class StepCountType, class ...Args>
  void discreteTimeJacobian(StepCountType currentStepNumber,
			    scalar_type time_np1,
			    scalar_type dt,
			    discrete_time_jacobian_type & J,
			    Args && ...args) const
  {
    preconditionable::discreteTimeJacobian(currentStepNumber, time_np1, dt,
					   J, std::forward<Args>(args)...);
    const auto & ynp1 = fomStatesMngr_(::pressio::ode::nPlusOne());
    precFunctor_(ynp1, time_np1, J);
  }

private:
  using preconditionable::fomStatesMngr_;

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  const PreconditionerType precFunctor_;
#else
  std::reference_wrapper<const PreconditionerType> precFunctor_;
#endif
};

template <class MaskerType, class maskable>
struct MaskDecoratorDiscreteTimeSystem : maskable
{
  using scalar_type = typename maskable::scalar_type;
  using state_type  = typename maskable::state_type;
  using discrete_time_residual_type = typename maskable::discrete_time_residual_type;
  using discrete_time_jacobian_type = typename maskable::discrete_time_jacobian_type;

  MaskDecoratorDiscreteTimeSystem() = delete;
  MaskDecoratorDiscreteTimeSystem(const MaskDecoratorDiscreteTimeSystem &) = default;
  MaskDecoratorDiscreteTimeSystem & operator=(const MaskDecoratorDiscreteTimeSystem &) = delete;
  MaskDecoratorDiscreteTimeSystem(MaskDecoratorDiscreteTimeSystem &&) = default;
  MaskDecoratorDiscreteTimeSystem & operator=(MaskDecoratorDiscreteTimeSystem &&) = delete;
  ~MaskDecoratorDiscreteTimeSystem() = default;

  template <class ... Args>
  MaskDecoratorDiscreteTimeSystem(
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
			pybind11::object maskerIn,
#else
			const MaskerType & maskerIn,
#endif
			Args && ... args)
    : maskable(std::forward<Args>(args)...),
      maskFunctor_(maskerIn),
      unmaskedResidual_(maskable::createDiscreteTimeResidual()),
      unmaskedJacobian_(maskable::createDiscreteTimeJacobian())
  {}

  discrete_time_residual_type createDiscreteTimeResidual() const{
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
    discrete_time_residual_type R(maskFunctor_.createApplyMaskResult(unmaskedResidual_));
#else
    discrete_time_residual_type R(maskFunctor_.get().createApplyMaskResult(unmaskedResidual_));
#endif
    return R;
  }

  discrete_time_jacobian_type createDiscreteTimeJacobian() const{
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
    discrete_time_jacobian_type J(maskFunctor_.createApplyMaskResult(unmaskedJacobian_));
#else
    discrete_time_jacobian_type J(maskFunctor_.get().createApplyMaskResult(unmaskedJacobian_));
#endif
    return J;
  }

  template<class StepCountType, class ...Args>
  void discreteTimeResidual(StepCountType currentStepNumber,
                            scalar_type time_np1,
                            scalar_type dt,
                            discrete_time_residual_type & R,
			    Args && ...args) const
  {
    maskable::discreteTimeResidual(currentStepNumber, time_np1, dt,
				   unmaskedResidual_, std::forward<Args>(args)...);
    maskFunctor_(unmaskedResidual_, time_np1, R);
  }

  template<class StepCountType, class ...Args>
  void discreteTimeJacobian(StepCountType currentStepNumber,
			    scalar_type time_np1,
			    scalar_type dt,
			    discrete_time_jacobian_type & J,
			    Args && ...args) const
  {
    maskable::discreteTimeJacobian(currentStepNumber, time_np1, dt,
				   unmaskedJacobian_, std::forward<Args>(args)...);
    maskFunctor_(unmaskedJacobian_, time_np1, J);
  }

private:
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  const MaskerType maskFunctor_;
#else
  std::reference_wrapper<const MaskerType> maskFunctor_;
#endif

  mutable discrete_time_residual_type unmaskedResidual_;
  mutable discrete_time_jacobian_type unmaskedJacobian_;
};

}}}}//end namespace pressio::rom::lspg::impl
#endif  // ROM_IMPL_ROM_LSPG_UNSTEADY_DISCRETE_TIME_DECORATORS_HPP_
