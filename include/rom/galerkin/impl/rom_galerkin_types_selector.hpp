/*
//@HEADER
// ************************************************************************
//
// rom_galerkin_types_selector.hpp
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

#ifndef ROM_GALERKIN_IMPL_ROM_GALERKIN_TYPES_SELECTOR_HPP_
#define ROM_GALERKIN_IMPL_ROM_GALERKIN_TYPES_SELECTOR_HPP_

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

template<typename T, class = void>
struct select_galerkin_types
{
  using residual_type = void;
  using jacobian_type = void;
};

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<typename T>
struct select_galerkin_types<
  T,
  mpl::enable_if_t<
    ::pressio::containers::predicates::is_dynamic_vector_wrapper_eigen<T>::value
    >
  >
{
  // for now use residual_type = state_type
  using residual_type = T;
  // the galerkin jacobian is a wrapper of eigen::matrix
  // for now make it column-major (which is default)
  using native_j_t = Eigen::MatrixXd;
  using jacobian_type = ::pressio::containers::DenseMatrix<native_j_t>;
};
#endif

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
template<typename T>
struct select_galerkin_types<
  T,
  mpl::enable_if_t<
    ::pressio::containers::predicates::is_rank1_tensor_wrapper_pybind<T>::value
    >
  >
{
  using scalar_t = typename ::pressio::containers::details::traits<T>::scalar_t;
  // for now use residual_type = state_type
  using residual_type = T;
  // the galerkin jacobian is a pybind11 tensor column-major
  using native_j_t = pybind11::array_t<scalar_t, pybind11::array::f_style>;
  using jacobian_type = ::pressio::containers::Tensor<2, native_j_t>;
};
#endif

}}}}
#endif  // ROM_GALERKIN_IMPL_ROM_GALERKIN_TYPES_SELECTOR_HPP_
