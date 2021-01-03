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

template <typename residual_type, typename fom_states_manager_t>
class ResidualPolicy
{
public:
  // typedef needed for decorators
  using data_type = residual_type;

public:
  ResidualPolicy() = delete;
  ResidualPolicy(const ResidualPolicy &) = default;
  ResidualPolicy & operator=(const ResidualPolicy &) = delete;
  ResidualPolicy(ResidualPolicy &&) = default;
  ResidualPolicy & operator=(ResidualPolicy &&) = delete;
  ~ResidualPolicy() = default;

  ResidualPolicy(fom_states_manager_t & fomStatesMngr)
    : fomStatesMngr_(fomStatesMngr){}

public:
  template <typename fom_system_t>
  residual_type create(const fom_system_t & fomSystemObj) const
  {
    return residual_type(fomSystemObj.createResidual());
  }

  template <typename lspg_state_t, typename fom_system_t>
  void compute(const lspg_state_t & romState,
	       residual_type & romResidual,
	       const fom_system_t & fomSystemObj) const
  {
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("lspg residual");
#endif

    fomStatesMngr_.get().reconstructCurrentFomState(romState);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->start("fom eval rhs");
#endif

    const auto & currFom = fomStatesMngr_.get().currentFomState();
    fomSystemObj.residual(*currFom.data(), *romResidual.data());

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("fom eval rhs");
    timer->stop("lspg residual");
#endif
  }

protected:
  std::reference_wrapper<fom_states_manager_t> fomStatesMngr_;
};

}}}}}
#endif  // ROM_LSPG_IMPL_STEADY_POLICIES_ROM_LSPG_STEADY_RESIDUAL_POLICY_HPP_
