/*
//@HEADER
// ************************************************************************
//
// rom_static_container_fom_states.hpp
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

#ifndef ROM_STATIC_CONTAINER_FOM_STATES_HPP_
#define ROM_STATIC_CONTAINER_FOM_STATES_HPP_

#include "rom_ConfigDefs.hpp"
#include "rom_fwd.hpp"
#include "../../containers/src/collection/containers_static_collection.hpp"
#include "../../ode/src/ode_states_container.hpp"

namespace pressio{ namespace rom{

template <typename fom_state_type, std::size_t n, typename reconstuctor_type>
class FomStatesStaticContainer<fom_state_type, n, reconstuctor_type>
{
  static_assert( ::pressio::containers::meta::is_wrapper<fom_state_type>::value
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
		 or containers::meta::is_array_pybind11<fom_state_type>::value
#endif
		 , "Currently, you can only create a FomStatesStaticContainer from types \
which have pressio wrappers.");

public:
  using data_type = ::pressio::containers::StaticCollection<fom_state_type, n>;

  FomStatesStaticContainer() = delete;

  template <typename ... Args>
  FomStatesStaticContainer(const reconstuctor_type & fomStateReconstr,
			   Args && ... args)
    : fomStateReconstrObj_(fomStateReconstr),
      data_( std::forward<Args>(args)... ){
    this->resetContainersToZero();
  }

  ~FomStatesStaticContainer() = default;

public:
  static constexpr std::size_t size() {
    return data_type::size();
  }

  const fom_state_type & getCRefToCurrentFomState() const{
    static_assert( n>=1, "Cannot call getCRefToCurrentFomState if n < 1");
    return data_[0];
  }

  const fom_state_type & getCRefToFomStatePrevStep() const{
    static_assert( n>=2, "Cannot call getCRefToFomStatePrevStep if n < 2");
    return data_[1];
  }

  const fom_state_type & getCRefToFomStatePrevPrevStep() const{
    static_assert( n>=3, "Cannot call getCRefToFomStatePrevPrevStep if n < 3");
    return data_[2];
  }


  template <typename rom_state_t>
  void reconstructCurrentFomState(const rom_state_t & romY)
  {
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("reconstruct fom state");
#endif

    static_assert( n>=1, "Cannot call reconstructCurrentFomState if n < 1");
    fomStateReconstrObj_(romY, data_[0]);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("reconstruct fom state");
#endif
  }


  template <std::size_t n2, typename rom_state_t>
  void reconstructFomOldStates(const ::pressio::ode::StatesContainer<rom_state_t, n2> & romYprev)
  {
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("reconstruct fom old state");
#endif

    static_assert( n>n2, "Cannot call reconstructFomOldStates if n > n2");
    for (std::size_t i=0; i<n2; i++){
      fomStateReconstrObj_(romYprev[i], data_[i+1]);
    }
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("reconstruct fom old state");
#endif
  }

private:
  /* set all entries to zero for all members */
  void resetContainersToZero(){
    for (std::size_t i=0; i<n; i++)
      ::pressio::containers::ops::set_zero(data_[i]);
  }

private:
  const reconstuctor_type & fomStateReconstrObj_  = {};

  // data[0] contains the current fom state, i.e. step = n
  // data[1] contains fom state at step n-1
  // data[2] contains fom state at n-2
  // etc..
  data_type data_;
};

}}//end namespace pressio::rom
#endif
