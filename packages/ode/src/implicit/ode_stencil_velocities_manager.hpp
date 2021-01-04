/*
//@HEADER
// ************************************************************************
//
// ode_stencil_velocities_manager.hpp
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

#ifndef ODE_IMPLICIT_ODE_STENCIL_VELOCITIES_MANAGER_HPP_
#define ODE_IMPLICIT_ODE_STENCIL_VELOCITIES_MANAGER_HPP_

namespace pressio{ namespace ode{ namespace implicitmethods{

namespace impl{
template<typename T, std::size_t n>
class StencilVelocitiesManager
{

public:
  using data_type = ::pressio::containers::IndexableStaticCollection<T, n>;

private:
  data_type rhs_;

public:
  StencilVelocitiesManager() = default;

  StencilVelocitiesManager(const T & oin)
    : rhs_(oin)
  {
    PRESSIOLOG_DEBUG("constructing from oin");
  }

  StencilVelocitiesManager(T && oin)
    : rhs_(std::move(oin))
  {
    PRESSIOLOG_DEBUG("constructing moving oin");
  }

  StencilVelocitiesManager(StencilVelocitiesManager const & other) = default;
  StencilVelocitiesManager & operator=(StencilVelocitiesManager const & other) = default;
  StencilVelocitiesManager(StencilVelocitiesManager && other) = default;
  StencilVelocitiesManager & operator=(StencilVelocitiesManager && other) = default;
  ~StencilVelocitiesManager() = default;

public:
  static constexpr std::size_t size(){ return n; }

  // for rhs we have np1, n, n-1, etc
  // because the velocities are handle by the stepper and owned inside
  // n+1
  template <std::size_t _n = n>
  mpl::enable_if_t<_n>=1, T & > velocityAt(ode::nPlusOne){ return rhs_(0); }

  template <std::size_t _n = n>
  mpl::enable_if_t<_n>=1, T const & > velocityAt(ode::nPlusOne) const { return rhs_(0); }

  template <std::size_t _n = n>
  mpl::enable_if_t<_n>=1, T & > operator()(ode::nPlusOne){ return rhs_(0); }

  template <std::size_t _n = n>
  mpl::enable_if_t<_n>=1, T const & > operator()(ode::nPlusOne) const { return rhs_(0); }

  // n
  template <std::size_t _n = n>
  mpl::enable_if_t<_n>=2, T & > velocityAt(ode::n){ return rhs_(1); }

  template <std::size_t _n = n>
  mpl::enable_if_t<_n>=2, T const & > velocityAt(ode::n) const { return rhs_(1); }

  template <std::size_t _n = n>
  mpl::enable_if_t<_n>=2, T & > operator()(ode::n){ return rhs_(1); }

  template <std::size_t _n = n>
  mpl::enable_if_t<_n>=2, T const & > operator()(ode::n) const { return rhs_(1); }

  // n-1
  template <std::size_t _n = n>
  mpl::enable_if_t<_n>=3, T & > velocityAt(ode::nMinusOne){ return rhs_(2); }

  template <std::size_t _n = n>
  mpl::enable_if_t<_n>=3, T const & > velocityAt(ode::nMinusOne) const { return rhs_(2); }

  template <std::size_t _n = n>
  mpl::enable_if_t<_n>=3, T & > operator()(ode::nMinusOne){ return rhs_(2); }

  template <std::size_t _n = n>
  mpl::enable_if_t<_n>=3, T const & > operator()(ode::nMinusOne) const { return rhs_(2); }

  // n-2
  template <std::size_t _n = n>
  mpl::enable_if_t<_n>=4, T & > velocityAt(ode::nMinusTwo){ return rhs_(3); }

  template <std::size_t _n = n>
  mpl::enable_if_t<_n>=4, T const & > velocityAt(ode::nMinusTwo) const { return rhs_(3); }

  template <std::size_t _n = n>
  mpl::enable_if_t<_n>=4, T & > operator()(ode::nMinusTwo){ return rhs_(3); }

  template <std::size_t _n = n>
  mpl::enable_if_t<_n>=4, T const & > operator()(ode::nMinusTwo) const { return rhs_(3); }
};
}// end namespace impl

template<typename T, std::size_t n>
using StencilVelocitiesManager = impl::StencilVelocitiesManager<T, n>;

}}}//end namespace pressio::ode::implicitmethods
#endif  // ODE_IMPLICIT_ODE_STENCIL_VELOCITIES_MANAGER_HPP_
