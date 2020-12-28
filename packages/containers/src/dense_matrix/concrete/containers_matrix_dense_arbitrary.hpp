/*
//@HEADER
// ************************************************************************
//
// containers_matrix_dense_arbitrary.hpp
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

#ifndef CONTAINERS_DENSE_MATRIX_CONCRETE_CONTAINERS_MATRIX_DENSE_ARBITRARY_HPP_
#define CONTAINERS_DENSE_MATRIX_CONCRETE_CONTAINERS_MATRIX_DENSE_ARBITRARY_HPP_

namespace pressio{ namespace containers{

template <typename wrapped_type>
class DenseMatrix<
  wrapped_type,
  mpl::enable_if_t<
    ::pressio::containers::predicates::is_admissible_as_dense_matrix_arbitrary<wrapped_type>::value
    >
  >
{
public:
  using this_t = DenseMatrix<wrapped_type>;
  using traits = details::traits<this_t>;
  using size_t = typename details::traits<this_t>::size_t;
  using sc_t   = typename details::traits<this_t>::scalar_t;

public:

  template<
    typename _wrapped_type = wrapped_type,
    mpl::enable_if_t<
      std::is_default_constructible<_wrapped_type>::value, int
    > = 0
  >
  DenseMatrix(){};

  template<
    typename _wrapped_type = wrapped_type,
    mpl::enable_if_t<
      std::is_constructible<_wrapped_type, size_t, size_t>::value, int
    > = 0
  >
  DenseMatrix(size_t nR, size_t nC) : data_(nR, nC){};

  explicit DenseMatrix(const wrapped_type & vecobj)
    : data_(vecobj){}

  DenseMatrix(DenseMatrix const & other)
    : data_(*other.data()){}

public:
  template<typename _wrapped_type = wrapped_type>
  mpl::enable_if_t<
  ::pressio::containers::predicates::has_method_extent<_wrapped_type>::value
  , size_t
  >
  extent(size_t k) const
  {
    assert( k==0 or k==1);
    return data_.extent(k);
  }

  template<typename _wrapped_type = wrapped_type>
  mpl::enable_if_t<
  ::pressio::containers::predicates::has_method_size_with_arg<_wrapped_type>::value
  , size_t
  >
  extent(size_t k) const{
    return data_.size(k);
  }

  sc_t & operator()(size_t i, size_t j){
    return data_(i, j);
  };

  sc_t const & operator()(size_t i, size_t j) const{
    return data_(i, j);
  };

  wrapped_type const * data() const{
    return &data_;
  }

  wrapped_type * data(){
    return &data_;
  }

private:
  wrapped_type data_ = {};

};//end class

}}//end namespace pressio::containers

#endif  // CONTAINERS_DENSE_MATRIX_CONCRETE_CONTAINERS_MATRIX_DENSE_ARBITRARY_HPP_
