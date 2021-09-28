/*
//@HEADER
// ************************************************************************
//
// rom_galerkin_projector_default.hpp
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

#ifndef ROM_IMPL_ROM_GALERKIN_PROJECTOR_DEFAULT_HPP_
#define ROM_IMPL_ROM_GALERKIN_PROJECTOR_DEFAULT_HPP_

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

/* default projector works for the non-hyperreduced case.
   It is supposed to compute "phi^T * operand" where :
   - operand = fom velocity (needed for explicit and implicit stepping)
     recall indeed that:
	dot(x_rom) = phi^T f(phi*x_rom)

   - operand = phi (needed for jacobian of residual for implicit stepping)
     recall indeed that for implicit time stepping we have:
	R = dot(x_rom) - phi^t f(phi*x_rom)
	so dR/dx_rom = ... - phi^T df/dx phi
	call B = df/dx phi

   so we need to compute below "phi^T f" and "phi^T B"
   Add overloads for phi and B of same type, if needed extend
*/

template <typename DecoderType>
struct DefaultProjector
{
  DefaultProjector() = delete;
  DefaultProjector(const DefaultProjector &) = default;
  DefaultProjector & operator=(const DefaultProjector &) = delete;
  DefaultProjector(DefaultProjector &&) = default;
  DefaultProjector & operator=(DefaultProjector &&) = delete;
  ~DefaultProjector() = default;

  DefaultProjector(const DecoderType & decoder)
    : decoderJacobian_(decoder.jacobianCRef()){}

  template<class OperandType, class TimeType, class ResultType>
  mpl::enable_if_t<::pressio::Traits<ResultType>::rank == 1>
  operator()(const OperandType & operand, const TimeType time, ResultType & result) const
  {
    (void)time;
    using scalar_t = typename ::pressio::Traits<ResultType>::scalar_type;
    using cnst = ::pressio::utils::Constants<scalar_t>;
    ::pressio::ops::product(::pressio::transpose(), cnst::one(),
			    decoderJacobian_.get(), operand,
			    cnst::zero(), result);
  }

  template<class OperandType, class TimeType, class ResultType>
  mpl::enable_if_t<::pressio::Traits<ResultType>::rank >= 2>
  operator()(const OperandType & operand, const TimeType time, ResultType & result) const
  {
    (void)time;
    using scalar_t = typename ::pressio::Traits<ResultType>::scalar_type;
    using cnst = ::pressio::utils::Constants<scalar_t>;
    ::pressio::ops::product(::pressio::transpose(), ::pressio::nontranspose(),
			    cnst::one(), decoderJacobian_.get(), operand,
			    cnst::zero(), result);
  }

private:
  std::reference_wrapper<const typename DecoderType::jacobian_type> decoderJacobian_;
};

}}}}//end  namespace pressio::rom::galerkin::impl
#endif  // ROM_IMPL_ROM_GALERKIN_PROJECTOR_DEFAULT_HPP_
