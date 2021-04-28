/*
//@HEADER
// ************************************************************************
//
// qr_out_of_place_base.hpp
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

#ifndef QR_BASE_QR_OUT_OF_PLACE_BASE_HPP_
#define QR_BASE_QR_OUT_OF_PLACE_BASE_HPP_

namespace pressio{ namespace qr{


template<typename derived_t, typename matrix_t, typename Q_t>
class QROutOfPlaceBase
  : private utils::details::CrtpBase<
  QROutOfPlaceBase<derived_t, matrix_t, Q_t>>
{

  using this_t = QROutOfPlaceBase<derived_t, matrix_t, Q_t>;

  /* workaround for nvcc issue with templates, see https://devtalk.nvidia.com/default/topic/1037721/nvcc-compilation-error-with-template-parameter-as-a-friend-within-a-namespace/ */
  template<typename DummyType> struct dummy{using type = DummyType;};
  friend typename dummy<derived_t>::type;

  friend utils::details::CrtpBase<this_t>;

public:
  void computeThin(const matrix_t & A){
    this->underlying().computeThinImpl(A);
  }

  const Q_t & cRefQFactor() const {
    return this->underlying().cRefQFactorImpl();
  }

  template <typename vec_in_t, typename vec_out_t>
  ::pressio::mpl::enable_if_t<
    containers::predicates::is_vector_wrapper<vec_in_t>::value and
    containers::predicates::is_vector_wrapper<vec_out_t>::value and
    meta::is_legitimate_vector_type_for_qr_project<vec_in_t, Q_t>::value
  >
  applyQTranspose(const vec_in_t & vecIn, vec_out_t & vecOut) const{
    this->underlying().applyQTransposeImpl(vecIn, vecOut);
  }

  template <typename vec_in_t, typename vec_out_t>
  ::pressio::mpl::enable_if_t<
    containers::predicates::is_vector_wrapper<vec_in_t>::value and
    containers::predicates::is_vector_wrapper<vec_out_t>::value
  >
  applyRTranspose(const vec_in_t & vecIn, vec_out_t & vecOut) const{
    this->underlying().applyRTransposeImpl(vecIn, vecOut);
  }

private:
  QROutOfPlaceBase() = default;
  ~QROutOfPlaceBase() = default;

};

}}//end namespace pressio::qr
#endif  // QR_BASE_QR_OUT_OF_PLACE_BASE_HPP_
