/*
//@HEADER
// ************************************************************************
//
// containers_multi_vector_view_column_vector_expressions.hpp
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

#ifndef CONTAINERS_MULTI_VECTOR_VIEW_COLUMN_VECTOR_EXPRESSIONS_HPP_
#define CONTAINERS_MULTI_VECTOR_VIEW_COLUMN_VECTOR_EXPRESSIONS_HPP_

namespace pressio{ namespace containers{ namespace expressions{

template <typename mv_t>
struct ViewColumnVectorExpr<
  mv_t, 
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_multi_vector_wrapper_eigen<mv_t>::value
    >
  >
  : public BaseExpr< ViewColumnVectorExpr<mv_t> >
{
  using native_t = typename ::pressio::containers::details::traits<mv_t>::wrapped_t;
  using scalar_t = typename ::pressio::containers::details::traits<mv_t>::scalar_t;

private:
  mv_t & mvObj_;
  const std::size_t vecIndex_;

public:
  ViewColumnVectorExpr() = delete;
  ~ViewColumnVectorExpr() = default;
  ViewColumnVectorExpr(const ViewColumnVectorExpr & other) = default;
  ViewColumnVectorExpr(ViewColumnVectorExpr && other) = default;
  ViewColumnVectorExpr & operator=(const ViewColumnVectorExpr & other) = default;
  ViewColumnVectorExpr & operator=(ViewColumnVectorExpr && other) = default;

  ViewColumnVectorExpr(mv_t & mvObjIn, const std::size_t vecIndexIn)
    : mvObj_(mvObjIn), vecIndex_(vecIndexIn){
    assert( vecIndexIn < mvObj_.numVectors() );
  }

  // this is for a column vector, so return the # of rows
  std::size_t extent(std::size_t i) const{ 
    assert(i==0);
    return mvObj_.extent(0);
  }

  const scalar_t & operator[](const std::size_t & rowIndex) const{
    return mvObj_(rowIndex, vecIndex_);
  }

  scalar_t & operator[](const std::size_t & rowIndex){
    return mvObj_(rowIndex, vecIndex_);
  }

  auto operator()()
    -> decltype( mvObj_.data()->col(vecIndex_) )
  {
    return mvObj_.data()->col(vecIndex_);
  }

  auto operator()() const
    -> decltype( std::declval<const native_t>().col(vecIndex_) )
  {
    return static_cast<native_t const *>(mvObj_.data())->col(vecIndex_);
  }
};


#ifdef PRESSIO_ENABLE_TPL_KOKKOS
template <typename mv_t>
struct ViewColumnVectorExpr<
  mv_t, 
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_multi_vector_wrapper_kokkos<mv_t>::value
    >
  >
  : public BaseExpr< ViewColumnVectorExpr<mv_t> >
{
private:
  const mv_t & mvObj_;
  const std::size_t vecIndex_;

public:
  ViewColumnVectorExpr() = delete;
  ~ViewColumnVectorExpr() = default;
  ViewColumnVectorExpr(const ViewColumnVectorExpr & other) = default;
  ViewColumnVectorExpr(ViewColumnVectorExpr && other) = default;
  ViewColumnVectorExpr & operator=(const ViewColumnVectorExpr & other) = default;
  ViewColumnVectorExpr & operator=(ViewColumnVectorExpr && other) = default;

  ViewColumnVectorExpr(const mv_t & mvObjIn, const std::size_t vecIndexIn)
    : mvObj_(mvObjIn), vecIndex_(vecIndexIn){}

  // this is for a column vector, so return the # of rows
  std::size_t extent(std::size_t i) const{ 
    assert(i==0);
    return mvObj_.extent(0);
  }

  auto operator()() const -> decltype( Kokkos::subview(*mvObj_.data(), Kokkos::ALL(), vecIndex_) )
  {
    return Kokkos::subview(*mvObj_.data(), Kokkos::ALL(), vecIndex_);
  }
};
#endif

}}} //end namespace pressio::containers::expressions

#endif