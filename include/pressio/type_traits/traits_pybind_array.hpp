/*
//@HEADER
// ************************************************************************
//
// traits_pybind_array.hpp
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

#ifndef TYPE_TRAITS_TRAITS_PYBIND_ARRAY_HPP_
#define TYPE_TRAITS_TRAITS_PYBIND_ARRAY_HPP_

namespace pressio{

template<class ScalarType>
struct Traits<pybind11::array_t<ScalarType, pybind11::array::c_style>>
  : public ::pressio::impl::ContainerTraits<
      -1, // no fixed rank
      ScalarType,
      std::size_t
    >
{

  static constexpr bool is_static = false;
  static constexpr bool is_dynamic = !is_static;

  static constexpr TensorIdentifier tensor_identifier
  = TensorIdentifier::Pybind;

  static constexpr int layout = 0;
  using scalar_type = ScalarType;
  using size_type = std::size_t;
};

template<class ScalarType>
struct Traits<pybind11::array_t<ScalarType, pybind11::array::f_style>>
  : public impl::ContainersSharedTraits<PackageIdentifier::Pybind, true, -1>
{

  static constexpr bool is_static = false;
  static constexpr bool is_dynamic = !is_static;

  static constexpr TensorIdentifier tensor_identifier
  = TensorIdentifier::Pybind;

  static constexpr int layout = 1;
  using scalar_type = ScalarType;
  using size_type = std::size_t;
};

}
#endif  // TYPE_TRAITS_TRAITS_PYBIND_ARRAY_HPP_
