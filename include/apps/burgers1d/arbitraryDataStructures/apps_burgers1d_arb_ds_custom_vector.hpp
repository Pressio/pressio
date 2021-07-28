/*
//@HEADER
// ************************************************************************
//
// apps_burgers1d_arb_ds_custom_vector.hpp
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

#ifndef APPS_BURGERS1D_ARBITRARYDATASTRUCTURES_APPS_BURGERS1D_ARB_DS_CUSTOM_VECTOR_HPP_
#define APPS_BURGERS1D_ARBITRARYDATASTRUCTURES_APPS_BURGERS1D_ARB_DS_CUSTOM_VECTOR_HPP_

namespace pressio{ namespace apps{ namespace arbds{

template <typename scalar_t>
class Vector
{
public:
  using value_type = scalar_t;
  using size_type  = std::size_t;
  using index_type = size_type;
  using data_type  = std::vector<scalar_t>;

public:
  Vector() = default;

  explicit Vector(std::size_t extent)
    : data_(extent){}

  void resize(size_type newSize){
    data_.resize(newSize);
  }

  size_type extent(size_type k) const{
    assert(k == 0);
    return data_.size();
  }

  value_type & operator()(size_type i){
    return data_[i];
  }
  value_type const & operator()(size_type i) const{
    return data_[i];
  }

  value_type & operator[](size_type i){
    return data_[i];
  }
  value_type const & operator[](size_type i) const{
    return data_[i];
  }

  data_type * data(){
    return &data_;
  }
  data_type const * data() const{
    return &data_;
  }

private:
  data_type data_;

};//end class

}}} //namespace pressio::apps::arbds
#endif  // APPS_BURGERS1D_ARBITRARYDATASTRUCTURES_APPS_BURGERS1D_ARB_DS_CUSTOM_VECTOR_HPP_
