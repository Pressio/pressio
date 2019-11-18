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

  static_assert( n>=1,
		 "You are trying to instantiate a state FomStatesContainer object with zero size.\
Something is not right, because this object would then be unusable.");

public:
  using data_type  = ::pressio::containers::StaticCollection<fom_state_type, n>;
  using value_type = fom_state_type;

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


  /* this method reconstructs the current FOM state */
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


  // we are overloading the left shift operator such that
  // it operates on a rom state object, reconstructs it and
  // updates the previous FOM states

  /* when n == 2, it means I only have current state and previous one
   * so when I need to reconstruct the previous state, I can simply
   * overwrite the data in data_[1] */
  template <
    typename rom_state_t,
    std::size_t _n = n,
    ::pressio::mpl::enable_if_t<_n == 2> * = nullptr
    >
  void operator << (const rom_state_t & romStateIn)
  {
    // then, reconstrct the FOM state at n-1
    fomStateReconstrObj_(romStateIn, data_[1]);
  }

  /* when n >= 3, we need to deep copy data to
   * such that y_t-2 goes into y_t-3,
   * and y_t-1 goes into y_t-2, etc.
   * and then finally we overwrite data_[0] */
  template <
    typename rom_state_t,
    std::size_t _n = n,
    ::pressio::mpl::enable_if_t<_n >= 3> * = nullptr
    >
  void operator << (const rom_state_t & romStateIn)
  {
    // copy all states back
    for (std::size_t i=n-2; i>=1; --i){
      const auto & src  = data_[i];
      auto & dest = data_[i+1];
      ::pressio::containers::ops::deep_copy(src, dest);
    }
    // then, reconstrct the FOM state at n-1
    fomStateReconstrObj_(romStateIn, data_[1]);
  }


//   template <std::size_t n2, typename rom_state_t>
//   void reconstructFomOldStates(const ::pressio::ode::StatesContainer<rom_state_t, n2> & romYprev)
//   {
// #ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
//     auto timer = Teuchos::TimeMonitor::getStackedTimer();
//     timer->start("reconstruct fom old state");
// #endif

//     static_assert( n>n2, "Cannot call reconstructFomOldStates if n <= n2");
//     for (std::size_t i=0; i<n2; i++){
//       fomStateReconstrObj_(romYprev[i], data_[i+1]);
//     }
// #ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
//     timer->stop("reconstruct fom old state");
// #endif
//   }

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
