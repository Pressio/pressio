/*
//@HEADER
// ************************************************************************
//
// rom_fwd.hpp
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

#ifndef ROM_FORWARD_DECLARATIONS_HPP_
#define ROM_FORWARD_DECLARATIONS_HPP_

#include "../../ode/src/ode_enum.hpp"

namespace pressio{ namespace rom{

template <
  typename fom_state_type, std::size_t n, typename reconstuctor_type,
  typename enable = void>
class FomStatesStaticContainer;

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
template <typename matrix_type, typename ops_t, typename enable = void>
struct PyLinearDecoder;
#endif

/* operators */
template<typename wrapped_type, typename enable = void>
class MultiVectorOperator;

template<typename wrapped_type, typename enable = void>
class MatrixOperator;


/*--- decorators ---*/
namespace decorator{

template <typename preconditionable, typename enable = void>
class Preconditioned;

template <typename maskable, typename enable = void>
class Masked;

}// namespace pressio::rom::decorator


namespace policy{

template <bool is_steady_problem>
struct QueryFomVelocityDefault;

template <bool is_steady_problem>
struct QueryFomApplyJacobianDefault;

struct QueryFomTimeDiscreteResidual;

}// namespace pressio::rom::policy


/* --- galerkin --- */
namespace galerkin{

template <
  typename fom_states_data_t,
  typename fom_rhs_t,
  typename decoder_jac_t,
  typename ud_ops = void
  >
class DefaultExplicitVelocityPolicy;

template <typename type_generator_t, typename enable = void>
class ProblemGenerator;

} // end namespace pressio::rom::galerkin


/*--- steady LSPG --- */
namespace lspg{ namespace steady{

template <
  typename fom_type, typename decoder_type, typename lspg_state_type,
  typename enable = void
  >
struct CommonTypes;

template <typename residual_type, typename fom_states_data_type, typename fom_rhs_eval_policy>
class ResidualPolicy;

template<
  typename fom_states_data, typename apply_jac_return_type,
  typename fom_apply_jac_policy, typename decoder_t
  >
class JacobianPolicy;

template <typename fom_type, typename decoder_type, typename lspg_state_type>
struct DefaultProblemType;

template <typename fom_type, typename decoder_type, typename lspg_state_type>
struct PreconditionedProblemType;

template<
  typename app_type,
  typename lspg_state_type,
  typename lspg_residual_type,
  typename lspg_jacobian_type,
  typename residual_policy_type,
  typename jacobian_policy_type,
  typename enable = void
  >
class System;

template <typename type_generator_t, typename enable = void>
class ProblemGenerator;

}} // end namespace lspg::steady


// These are here for backward compatibility, should be deleted at some point
template <typename ... Args>
using DefaultLSPGSteadyTypeGenerator = ::pressio::rom::lspg::steady::DefaultProblemType<Args...>;

template <typename ... Args>
using PreconditionedLSPGSteadyTypeGenerator = ::pressio::rom::lspg::steady::PreconditionedProblemType<Args...>;

template <typename problem_type, typename enable = void>
using LSPGSteadyProblemGenerator = ::pressio::rom::lspg::steady::ProblemGenerator<problem_type>;


}} // end namespace pressio::rom
#endif
