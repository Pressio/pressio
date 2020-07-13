/*
//@HEADER
// ************************************************************************
//
// containers_matrix_dense_distributed_epetra.hpp
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

#ifndef CONTAINERS_MATRIX_CONCRETE_CONTAINERS_MATRIX_DENSE_DISTRIBUTED_EPETRA_HPP_
#define CONTAINERS_MATRIX_CONCRETE_CONTAINERS_MATRIX_DENSE_DISTRIBUTED_EPETRA_HPP_

namespace pressio{ namespace containers{

template <typename wrapped_type>
class Matrix<wrapped_type,
	     ::pressio::mpl::enable_if_t<
	       containers::predicates::is_dense_matrix_epetra<
		 wrapped_type>::value>
	     >
  : public MatrixDistributedBase< Matrix<wrapped_type> >
{

  using this_t = Matrix<wrapped_type>;
  using traits_t = details::traits<this_t>;
  using sc_t = typename traits_t::scalar_t;
  using LO_t = typename traits_t::local_ordinal_t;
  using GO_t = typename traits_t::global_ordinal_t;
  using comm_t =  typename traits_t::communicator_t;
  using wrap_t = typename traits_t::wrapped_t;
  using row_map_t = typename traits_t::row_map_t;

public:

  Matrix() = delete;

  Matrix(const row_map_t & rowMap, GO_t ncols)
    : data_(rowMap, ncols){}

  explicit Matrix(const wrap_t & objin)
    : data_(objin){}

  ~Matrix() = default;

public:
  sc_t & operator()(LO_t irow, GO_t icol){
    assert(icol < this->globalCols() );
    assert(irow < this->localRows() );
    return data_[icol][irow];
  }

  sc_t const & operator()(LO_t irow, GO_t icol)const{
    assert(icol < this->globalCols() );
    assert(irow < this->localRows() );
    return data_[icol][irow];
  }

public:
  wrap_t const * data() const{
    return &data_;
  }

  wrap_t * data(){
    return &data_;
  }

private:
  friend MatrixDistributedBase< this_t >;
  wrap_t data_ = {};

};//end class

}}//end namespace pressio::containers
#endif  // CONTAINERS_MATRIX_CONCRETE_CONTAINERS_MATRIX_DENSE_DISTRIBUTED_EPETRA_HPP_
