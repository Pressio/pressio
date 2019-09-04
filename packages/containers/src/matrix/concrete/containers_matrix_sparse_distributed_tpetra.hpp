/*
//@HEADER
// ************************************************************************
//
// containers_matrix_sparse_distributed_tpetra.hpp
//                     		      Pressio 
// Copyright 2019 National Technology & Engineering Solutions of Sandia,LLC 
//							      (NTESS)
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

#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_MATRIX_CONCRETE_MATRIX_SPARSE_DISTRIBUTED_TPETRA_HPP_
#define CONTAINERS_MATRIX_CONCRETE_MATRIX_SPARSE_DISTRIBUTED_TPETRA_HPP_

#include "../../shared_base/containers_container_base.hpp"
#include "../../shared_base/containers_container_distributed_mpi_base.hpp"
#include "../../shared_base/containers_container_distributed_trilinos_base.hpp"

#include "../base/containers_matrix_base.hpp"
#include "../base/containers_matrix_sparse_base.hpp"
#include "../base/containers_matrix_distributed_base.hpp"
#include "../base/containers_matrix_sparse_distributed_base.hpp"
#include "../base/containers_matrix_sparse_distributed_trilinos_base.hpp"

namespace pressio{ namespace containers{

template <typename wrapped_type>
class Matrix<wrapped_type,
	     ::pressio::mpl::enable_if_t<
	       containers::meta::is_sparse_matrix_tpetra<
		 wrapped_type>::value>
	     >
  : public ContainerBase< Matrix<wrapped_type>, wrapped_type >,
    public MatrixDistributedBase< Matrix<wrapped_type> >,
    public MatrixSparseDistributedTrilinosBase< Matrix<wrapped_type> >,
    public ContainerDistributedMpiBase< Matrix<wrapped_type>,
      typename details::traits<Matrix<wrapped_type>>::communicator_t>/*
    public MatrixBase< Matrix<wrapped_type> >,
    public MatrixSparseBase< Matrix<wrapped_type> >,
    public MatrixSparseDistributedBase< Matrix<wrapped_type> >,
    */
{

  using derived_t = Matrix<wrapped_type>;
  using traits_t = details::traits<derived_t>;
  using sc_t = typename traits_t::scalar_t;
  using LO_t = typename traits_t::local_ordinal_t;
  using GO_t = typename traits_t::global_ordinal_t;
  using comm_t =  typename traits_t::communicator_t;
  using wrap_t = typename traits_t::wrapped_t;
  using row_map_t = typename traits_t::row_map_t;
  using col_map_t = typename traits_t::col_map_t;
  using range_map_t = typename traits_t::range_map_t;
  using domain_map_t = typename traits_t::domain_map_t;

public:
  Matrix() = delete;

  Matrix(Teuchos::RCP<const row_map_t>  rowMap,
	 LO_t NumEntriesPerRow)
    : data_(rowMap, NumEntriesPerRow){}

  Matrix(const derived_t &) = delete;

  // explicit Matrix(const wrap_t & objin)
  //   : data_(objin.getRowMap(),
  // 	    objin.getColMap(),
  // 	    objin.getGlobalMaxNumRowEntries()){
  //   assert(objin.isFillComplete());
  // }

  ~Matrix() = default;

private:
  wrap_t const * dataImpl() const{
    return &data_;
  }

  wrap_t * dataImpl(){
    return &data_;
  }

  wrap_t dataCpImpl(){
    return data_;
  }

  LO_t localRowsImpl() const{
    return data_.getNodeNumRows();
  }

  LO_t localColsImpl() const{
    return data_.getNodeNumCols();
  }

  GO_t globalRowsImpl() const{
    return data_.getGlobalNumRows();
  }

  GO_t globalColsImpl() const{
    return data_.getGlobalNumCols();
  }

  comm_t commImpl() const{
    return data_.getComm();
  }

  bool isFillingCompletedImpl() const{
    return data_.isFillComplete();
  }

  void fillingIsCompletedImpl(){
    data_.fillComplete();
  }

  void fillingIsCompletedImpl(Teuchos::RCP<const domain_map_t> const & dmap,
                              Teuchos::RCP<const range_map_t> const & rmap){
    // this is needed for rectangular matrices
    data_.fillComplete(dmap, rmap);
  }

private:
  friend ContainerBase< derived_t, wrapped_type >;
  friend MatrixDistributedBase< derived_t >;
  friend MatrixSparseDistributedTrilinosBase< derived_t >;
  friend ContainerDistributedMpiBase< derived_t, comm_t >;
  // friend MatrixBase< derived_t >;
  // friend MatrixSparseBase< derived_t >;
  // friend MatrixSparseDistributedBase< derived_t >;

private:
  wrap_t data_ = {};

};//end class

}}//end namespace pressio::containers
#endif
#endif
