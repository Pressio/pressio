/*
//@HEADER
// ************************************************************************
//
// containers_matrix_sparse_distributed_trilinos_base.hpp
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

#ifndef CONTAINERS_MATRIX_BASE_MATRIX_SPARSE_DISTRIBUTED_TRILINOS_BASE_HPP_
#define CONTAINERS_MATRIX_BASE_MATRIX_SPARSE_DISTRIBUTED_TRILINOS_BASE_HPP_

#include "../containers_matrix_traits.hpp"

namespace pressio{
namespace containers{

template<typename derived_type>
class MatrixSparseDistributedTrilinosBase
  : private utils::details::CrtpBase<
  MatrixSparseDistributedTrilinosBase<derived_type>>{

  static_assert( details::traits<derived_type>::is_shared_mem==0,
  "OOPS: non-distributed matrix inheriting from sparse distributed trilinos!");

  using traits_t = details::traits<derived_type>;
  using row_map_t = typename traits_t::row_map_t;
  using col_map_t = typename traits_t::col_map_t;
  using range_map_t = typename traits_t::range_map_t;
  using domain_map_t = typename traits_t::domain_map_t;

public:

  bool isFillingCompleted() const{
    return this->underlying().isFillingCompletedImpl();}

  void fillingIsCompleted(){
    this->underlying().fillingIsCompletedImpl();}

  void fillingIsCompleted(domain_map_t const & dmap,
			  range_map_t const & rmap){
    this->underlying().fillingIsCompletedImpl(dmap, rmap);}

  range_map_t const & getRangeDataMap() const{
    assert(this->isFillingCompleted());
    return this->underlying().getRangeDataMapImpl();}

  domain_map_t const & getDomainDataMap() const{
    assert(this->isFillingCompleted());
    return this->underlying().getDomainDataMapImpl();}

  row_map_t const & getRowDataMap() const{
    return this->underlying().getRowDataMapImpl();}

  col_map_t const & getColDataMap() const{
    return this->underlying().getColDataMapImpl();}

  bool hasSameRangeDataMapAs(derived_type const & other) const{
    return this->underlying().hasSameRangeDataMapAsImpl(other);}

  bool hasSameDomainDataMapAs(derived_type const & other) const{
    return this->underlying().hasSameDomainDataMapAsImpl(other);}

  bool hasSameRowDataMapAs(derived_type const & other) const{
    return this->underlying().hasSameRowDataMapAsImpl(other);}

  bool hasSameColDataMapAs(derived_type const & other) const{
    return this->underlying().hasSameColDataMapAsImpl(other);}

private:
  /* workaround for nvcc issue with templates, see https://devtalk.nvidia.com/default/topic/1037721/nvcc-compilation-error-with-template-parameter-as-a-friend-within-a-namespace/ */
  template<typename DummyType> struct dummy{using type = DummyType;};
  friend typename dummy<derived_type>::type;

  friend utils::details::CrtpBase<
    MatrixSparseDistributedTrilinosBase<derived_type>>;

  MatrixSparseDistributedTrilinosBase() = default;
  ~MatrixSparseDistributedTrilinosBase() = default;

};//end class

} // end namespace containers
}//end namespace pressio
#endif
