/*
//@HEADER
// ************************************************************************
//
// ode_aux_states_manager.hpp
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

#ifndef ODE_ODE_AUX_STATES_MANAGER_HPP_
#define ODE_ODE_AUX_STATES_MANAGER_HPP_

namespace pressio{ namespace ode{

template<typename T, std::size_t n>
class AuxStatesManager;

template<typename T, std::size_t n>
class AuxStatesManager
{
  static_assert
  (::pressio::containers::predicates::is_wrapper<T>::value,
   "AuxStatesManager only supports pressio containers.");

  static_assert
  (!::pressio::containers::predicates::is_expression<T>::value,
   "AuxStatesManager does NOT support pressio expressions.");

public:
  using data_type = ::pressio::containers::IndexableStaticCollection<T, n>;

  template <
    typename _T = T,
    mpl::enable_if_t<
      std::is_default_constructible<_T>::value, int
      > = 0
  >
  AuxStatesManager(){};

  template <
    typename ... Args,
    mpl::enable_if_t<sizeof...(Args) >= 1, int > = 0
    >
  AuxStatesManager(Args && ... args)
    : data_( std::forward<Args>(args)... ){}

  // copy cnstr
  AuxStatesManager(AuxStatesManager const & other) = default;
  // copy assignment
  AuxStatesManager & operator=(AuxStatesManager const & other) = default;
  // move cnstr
  AuxStatesManager(AuxStatesManager && other) = default;
  // move assignment
  AuxStatesManager & operator=(AuxStatesManager && other) = default;

  ~AuxStatesManager() = default;

public:
  static constexpr std::size_t size() {
    return data_type::size();
  }

  // n-1
  T & refStateAt(ode::nMinusOne tag){ return data_(0); }
  T const & cRefStateAt(ode::nMinusOne tag) const{ return data_(0); }
  T & stateAt(ode::nMinusOne tag){ return data_(0); }
  T const & stateAt(ode::nMinusOne tag) const{ return data_(0); }

  // n-2
  T & refStateAt(ode::nMinusTwo tag){ return data_(1); }
  T const & cRefStateAt(ode::nMinusTwo tag) const{ return data_(1); }
  T & stateAt(ode::nMinusTwo tag){ return data_(1); }
  T const & stateAt(ode::nMinusTwo tag) const{ return data_(1); }

  // n-3
  T & refStateAt(ode::nMinusThree tag){ return data_(2); }
  T const & cRefStateAt(ode::nMinusThree tag) const{ return data_(2); }
  T & stateAt(ode::nMinusThree tag){ return data_(2); }
  T const & stateAt(ode::nMinusThree tag) const{ return data_(2); }

  // n-4
  T & refStateAt(ode::nMinusFour tag){ return data_(3); }
  T const & cRefStateAt(ode::nMinusFour tag) const{ return data_(3); }
  T & stateAt(ode::nMinusFour tag){ return data_(3); }
  T const & stateAt(ode::nMinusFour tag) const{ return data_(3); }

private:
  data_type data_;
};

}}//end namespace pressio::ode
#endif  // ODE_ODE_AUX_STATES_MANAGER_HPP_
