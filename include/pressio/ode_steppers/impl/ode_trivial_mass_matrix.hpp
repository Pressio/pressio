/*
//@HEADER
// ************************************************************************
//
// ode_explicit_stepper.hpp
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

#ifndef ODE_STEPPERS_IMPL_ODE_TRIVIAL_MASS_MATRIX_HPP_
#define ODE_STEPPERS_IMPL_ODE_TRIVIAL_MASS_MATRIX_HPP_

namespace pressio{ namespace ode{ namespace impl{

struct NoOpMassMatrix{};

template<class T> struct is_trivial_mass_matrix : std::false_type{};
template<> struct is_trivial_mass_matrix<NoOpMassMatrix> : std::true_type{};

template<class T>
struct createMassMatrix{
  template<class SysType>
  T operator()(const SysType & system){
    return system.createMassMatrix();
  }
};

template<>
struct createMassMatrix<NoOpMassMatrix>{
  template<class SysType>
  NoOpMassMatrix operator()(const SysType & /*system*/){
    return NoOpMassMatrix();
  }
};

template<class T, class = void>
struct find_mass_matrix_if_any_or_noop{
  using type = NoOpMassMatrix;
};

template<class T>
struct find_mass_matrix_if_any_or_noop<
  T, mpl::enable_if_t<system_has_complete_mass_matrix_api<T>::value>
  >{
  using type = typename T::mass_matrix_type;
};

}}}//end namespace pressio::ode::explicitmethods::impl
#endif  // ODE_STEPPERS_IMPL_ODE_EXPLICIT_STEPPER_HPP_
