/*
//@HEADER
// ************************************************************************
//
// containers_matrix_sparse_distributed_epetra.hpp
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
#ifndef CONTAINERS_MATRIX_CONCRETE_MATRIX_SPARSE_DISTRIBUTED_EPETRA_HPP_
#define CONTAINERS_MATRIX_CONCRETE_MATRIX_SPARSE_DISTRIBUTED_EPETRA_HPP_

#include "../../shared_base/containers_container_base.hpp"
#include "../../shared_base/containers_container_distributed_mpi_base.hpp"

#include "../base/containers_matrix_base.hpp"
#include "../base/containers_matrix_sparse_base.hpp"
#include "../base/containers_matrix_distributed_base.hpp"
#include "../base/containers_matrix_sparse_distributed_base.hpp"
#include "../base/containers_matrix_sparse_distributed_trilinos_base.hpp"

namespace pressio{ namespace containers{

template <typename wrapped_type>
class Matrix<wrapped_type,
	     ::pressio::mpl::enable_if_t<
	       containers::meta::is_sparse_matrix_epetra<
		 wrapped_type>::value>
	     >
  : public ContainerBase< Matrix<wrapped_type>, wrapped_type >,
    public MatrixBase< Matrix<wrapped_type> >,
    public MatrixSparseBase< Matrix<wrapped_type> >,
    public MatrixDistributedBase< Matrix<wrapped_type> >,
    public MatrixSparseDistributedBase< Matrix<wrapped_type> >,
    public MatrixSparseDistributedTrilinosBase< Matrix<wrapped_type> >,
    public ContainerDistributedMpiBase< Matrix<wrapped_type>,
      typename details::traits<Matrix<wrapped_type>>::communicator_t >
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
  using crs_graph_t = typename traits_t::crs_graph_t;

public:
  Matrix() = delete;

  Matrix(const row_map_t & rowMap,
		  LO_t NumEntriesPerRow,
		  bool StaticProfile=false)
    : data_(Epetra_DataAccess::Copy, rowMap,
	    NumEntriesPerRow, StaticProfile){}

  // template <typename T,
  // 	    typename std::enable_if<
  // 	      meta::publicly_inherits_from<T, Epetra_BlockMap>::value
  // 	      >::type * = nullptr>
  // explicit Matrix(const T & rowMap,
  // 		  LO_t NumEntriesPerRow,
  // 		  bool StaticProfile=false)
  //   : data_(Epetra_DataAccess::Copy, static_cast<const Epetra_Map &>(rowMap),
  // 	    NumEntriesPerRow, StaticProfile){}

  explicit Matrix(const wrap_t & objin)
    : data_(objin){
    // cannot copy construct Crs matrix is filling is not completed
    assert(this->isFillingCompleted());
  }

  ~Matrix() = default;

private:
  wrap_t const * dataImpl() const{
    return &data_;
  }
  wrap_t * dataImpl(){
    return &data_;
  }

  void matchLayoutWithImpl(const derived_t & other){
    data_.ReplaceRowMap( other.getRowDataMap() );
    if ( other.data()->HaveColMap() )
      data_.ReplaceColMap( other.getColDataMap() );
  }

  LO_t localRowsImpl() const{
    return data_.NunMyRows();
  }

  LO_t localColsImpl() const{
    return data_.NunMyCols();
  }

  GO_t globalRowsImpl() const{
    return data_.NumGlobalRows();
  }

  GO_t globalColsImpl() const{
    return data_.NumGlobalCols();
  }

  comm_t const & commCRefImpl() const{
    return data_.Comm();
  }

  bool isDistributedGloballyImpl() const{
    return data_.DistributedGlobal();
  }

  bool isFillingCompletedImpl() const{
    return data_.Filled();
  }

  void fillingIsCompletedImpl(){
    data_.FillComplete();
  }

  void fillingIsCompletedImpl(domain_map_t const & dmap,
			      range_map_t const & rmap){
    // this is needed for rectangular matrices
    /* also note that epetra crs matrix ONLY
       takes Epetra_Map as map type.
       There is no way (that I can see)
       to pass a Epetra_BlockMap */
    data_.FillComplete(dmap, rmap);
  }

  row_map_t const & getRowDataMapImpl() const{
    return data_.RowMap();
  }

  col_map_t const & getColDataMapImpl() const{
    return data_.ColMap();
  }

  range_map_t const & getRangeDataMapImpl() const{
    return data_.RangeMap();
  }

  domain_map_t const & getDomainDataMapImpl() const{
    return data_.DomainMap();
  }

  bool hasSameRangeDataMapAsImpl(derived_t const & other) const{
    return data_.RangeMap().SameAs(other.getRangeDataMap());
  }

  bool hasSameDomainDataMapAsImpl(derived_t const & other) const{
    return data_.DomainMap().SameAs(other.getDomainDataMap());
  }

  bool hasSameRowDataMapAsImpl(derived_t const & other) const{
    return data_.RowMap().SameAs(other.getRowDataMap());
  }

  bool hasSameColDataMapAsImpl(derived_t const & other) const{
    return data_.ColMap().SameAs(other.getColDataMap());
  }

  void insertGlobalValuesImpl(GO_t targetRow,
			      LO_t numEntries,
			      const sc_t * values,
			      const GO_t * indices)
  {
    data_.InsertGlobalValues(targetRow,
			     numEntries,
			     values,
			     indices);
  }

  void addToDiagonalImpl(sc_t value) {
    const auto &  myRowMap = this->getRowDataMap();
    // get current diagonal
    Epetra_Vector diag(myRowMap);
    data_.ExtractDiagonalCopy(diag);
    // add to diagonal
    for (auto i=0; i<myRowMap.NumMyElements(); i++)
      diag[i] += value;
    // replace
    data_.ReplaceDiagonalValues(diag);
  }

  void scaleImpl(sc_t value) {
    data_.Scale(value);
  }

private:
  friend ContainerBase< derived_t, wrapped_type >;
  friend MatrixBase< derived_t >;
  friend MatrixSparseBase< derived_t >;
  friend MatrixDistributedBase< derived_t >;
  friend MatrixSparseDistributedBase< derived_t >;
  friend MatrixSparseDistributedTrilinosBase< derived_t >;
  friend ContainerDistributedMpiBase< derived_t, comm_t >;

private:
  wrap_t data_ = {};

};//end class

}}//end namespace pressio::containers
#endif
#endif
