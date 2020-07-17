/*
//@HEADER
// ************************************************************************
//
// rom_lspg_steady_residual_policy.hpp
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

#ifndef ROM_LSPG_IMPL_STEADY_POLICIES_ROM_LSPG_STEADY_RESIDUAL_POLICY_HPP_
#define ROM_LSPG_IMPL_STEADY_POLICIES_ROM_LSPG_STEADY_RESIDUAL_POLICY_HPP_

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{ namespace steady{

template <typename residual_type, typename fom_states_manager_t, typename ud_ops_type>
class ResidualPolicy
{
public:
  using residual_t = residual_type;
  using ud_ops_t = ud_ops_type;

public:
  ResidualPolicy() = delete;
  ~ResidualPolicy() = default;

  // 1. void ops
  template <
    typename _fom_states_manager_t = fom_states_manager_t,
    typename _ud_ops_t = ud_ops_type,
    ::pressio::mpl::enable_if_t< std::is_void<_ud_ops_t>::value, int > = 0
    >
  ResidualPolicy(_fom_states_manager_t & fomStatesMngr)
    : fomStatesMngr_(fomStatesMngr){}

  // 2. nonvoid ops
  template <
    typename _fom_states_manager_t = fom_states_manager_t,
    typename _ud_ops_t = ud_ops_type,
    ::pressio::mpl::enable_if_t< !std::is_void<_ud_ops_t>::value, int > = 0
    >
  ResidualPolicy(_fom_states_manager_t & fomStatesMngr, const _ud_ops_t & udOps)
    : fomStatesMngr_(fomStatesMngr), udOps_{&udOps}{}

public:
  template <typename fom_system_t>
  mpl::enable_if_t< !::pressio::ops::predicates::is_object_pybind<fom_system_t>::value, residual_t >
  create(const fom_system_t & fomSystemObj) const
  {
    return residual_t(fomSystemObj.createResidual());
  }

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  template <typename fom_system_t>
  mpl::enable_if_t< ::pressio::ops::predicates::is_object_pybind<fom_system_t>::value, residual_t >
  create(const fom_system_t & fomSystemObj) const
  {
    const auto & currentFom = fomStatesMngr_.getCRefToCurrentFomState();
    return residual_t(fomSystemObj.attr("residual")(*currentFom.data()));
  }
#endif

  template <typename lspg_state_t, typename fom_system_t, typename norm_value_type>
  void compute(const lspg_state_t & romState,
	       residual_type & romResidual,
	       const fom_system_t & fomSystemObj,
	       ::pressio::Norm normKind,
	       norm_value_type & normValue) const
  {
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("lspg residual");
#endif

    fomStatesMngr_.reconstructCurrentFomState(romState);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->start("fom eval rhs");
#endif

    ::pressio::rom::queryFomResidual(fomSystemObj, fomStatesMngr_.getCRefToCurrentFomState(), romResidual);

    if (normKind == ::pressio::Norm::L2)
      normValue = ::pressio::ops::norm2(romResidual);
    else if (normKind == ::pressio::Norm::L1)
      normValue = ::pressio::ops::norm1(romResidual);
    else
      throw std::runtime_error("Invalid norm kind for lspg unsteady residual policy");

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("fom eval rhs");
    timer->stop("lspg residual");
#endif
  }


protected:
  fom_states_manager_t & fomStatesMngr_;

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  // here we do this conditional type because it seems when ud_ops_t= pybind11::object
  // it only works if we copy the object. Need to figure out if we can leave ptr in all cases.
  typename std::conditional<
    ::pressio::mpl::is_same<ud_ops_t, pybind11::object>::value, ud_ops_t,
    const ud_ops_t *
    >::type udOps_ = {};
#else
  const ud_ops_t * udOps_ = {};
#endif

};//end class

}}}}}
#endif  // ROM_LSPG_IMPL_STEADY_POLICIES_ROM_LSPG_STEADY_RESIDUAL_POLICY_HPP_
