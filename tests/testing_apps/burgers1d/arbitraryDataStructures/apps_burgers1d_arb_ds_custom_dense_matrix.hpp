/*
//@HEADER
// ************************************************************************
//
// apps_burgers1d_arb_ds_custom_dense_matrix.hpp
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

#ifndef APPS_BURGERS1D_ARBITRARYDATASTRUCTURES_APPS_BURGERS1D_ARB_DS_CUSTOM_DENSE_MATRIX_HPP_
#define APPS_BURGERS1D_ARBITRARYDATASTRUCTURES_APPS_BURGERS1D_ARB_DS_CUSTOM_DENSE_MATRIX_HPP_

namespace pressio{ namespace apps{ namespace arbds{

template <typename scalar_t>
class DenseMatrix
{
public:
  using value_type = scalar_t;
  using size_type  = std::size_t;
  using index_type = size_type;

  using data_type  = std::vector<std::vector<scalar_t>>;
public:
  // DenseMatrix() = default;

  DenseMatrix(index_type nRows, index_type nCols)
    : data_(nRows), numRows_{nRows}, numCols_{nCols}{
      for (auto & it : data_) it.resize(nCols);
    }

  void resize(index_type newRows, index_type newCols){
    data_.resize(newRows);
    for (auto & it : data_) it.resize(newCols);
  }

  index_type extent(index_type k) const{
    assert(k <= 1);
    assert(data_.empty() == false);

    return (k==0) ? numRows_ : numCols_;
  }

  value_type & operator()(index_type i, index_type j){
    return data_[i][j];
  }
  value_type const & operator()(index_type i, index_type j) const{
    return data_[i][j];
  }

  data_type * data(){
    return &data_;
  }
  data_type const * data() const{
    return &data_;
  }

private:
  data_type data_;
  size_type numRows_ {};
  size_type numCols_ {};

};//end class

}}} //namespace pressio::apps::arbds
#endif  // APPS_BURGERS1D_ARBITRARYDATASTRUCTURES_APPS_BURGERS1D_ARB_DS_CUSTOM_DENSE_MATRIX_HPP_
