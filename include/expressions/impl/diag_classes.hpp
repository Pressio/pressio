/*
//@HEADER
// ************************************************************************
//
// containers_diag_classes.hpp
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

#ifndef CONTAINERS_EXPRESSIONS_DIAG_CONTAINERS_DIAG_CLASSES_HPP_
#define CONTAINERS_EXPRESSIONS_DIAG_CONTAINERS_DIAG_CLASSES_HPP_

namespace pressio{ namespace expressions{ namespace impl{

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template <typename MatrixType>
struct DiagExpr<
  MatrixType,
  ::pressio::mpl::enable_if_t<
    ::pressio::is_dense_matrix_eigen<MatrixType>::value
    >
  >
{
  using this_t = DiagExpr<MatrixType>;
  using mytraits = diag_traits<this_t>;
  using size_t = typename mytraits::size_type;
  using ref_t = typename mytraits::reference_type;
  using const_ref_t = typename mytraits::const_reference_type;
  using native_expr_t = typename mytraits::native_expr_type;

private:
  std::reference_wrapper<MatrixType> matObj_;
  native_expr_t nativeExprObj_;
  size_t numRows_ = {};
  size_t numCols_ = {};
  size_t extent_ = {};

public:
  DiagExpr() = delete;

  DiagExpr(const DiagExpr & other) = default;
  DiagExpr & operator=(const DiagExpr & other) = delete;

  DiagExpr(DiagExpr && other) = default;
  DiagExpr & operator=(DiagExpr && other) = delete;
  ~DiagExpr() = default;

  DiagExpr(MatrixType & matObjIn)
    : matObj_(matObjIn),
      nativeExprObj_(matObj_.get().diagonal()),
      numRows_(matObj_.get().rows()),
      numCols_(matObj_.get().cols()),
      extent_(matObj_.get().rows())
  {
    assert(numRows_ == numCols_);
  }

  size_t extent(size_t i) const{
    assert(i==0);
    return extent_;
  }

  native_expr_t const & native() const{
    return nativeExprObj_;
  }

  native_expr_t & native(){
    return nativeExprObj_;
  }

  ref_t operator()(size_t i)
  {
    assert(i < (size_t)extent_);
    return nativeExprObj_(i);
  }

  const_ref_t operator()(size_t i) const
  {
    assert(i < (size_t)extent_);
    return nativeExprObj_(i);
  }
};
#endif


#ifdef PRESSIO_ENABLE_TPL_KOKKOS
template <typename MatrixType>
struct DiagExpr<
  MatrixType,
  ::pressio::mpl::enable_if_t<
    ::pressio::is_dense_matrix_kokkos<MatrixType>::value
    >
  >
{
  using this_t		= DiagExpr<MatrixType>;
  using mytraits	= diag_traits<this_t>;
  using size_t		= typename mytraits::size_type;
  using ref_t		= typename mytraits::reference_type;
  using native_expr_t	= typename mytraits::native_expr_type;

private:
  std::reference_wrapper<MatrixType> matObj_;
  native_expr_t nativeExprObj_;
  size_t extent_ = {};

  using natexpr_layout = typename native_expr_t::traits::array_layout;
  // for now leave this assert, then remove later
  static_assert
  (std::is_same<natexpr_layout, Kokkos::LayoutStride>::value,
   "The layout for the native type of the diagonal kokkos expression does not match the strided layout expected");

public:
  DiagExpr() = delete;
  DiagExpr(const DiagExpr & other) = default;
  DiagExpr & operator=(const DiagExpr & other) = delete;
  DiagExpr(DiagExpr && other) = default;
  DiagExpr & operator=(DiagExpr && other) = delete;
  ~DiagExpr() = default;

  DiagExpr(MatrixType & M)
    : matObj_(M),
      nativeExprObj_(M.data(), natexpr_layout(M.extent(0), M.stride(0)+M.stride(1))),
      extent_(M.extent(0))
  {
    // make sure the diagonal is taken on a square matrix
    assert(M.extent(0) == M.extent(1));
  }

public:
  size_t extent(size_t i) const{
    assert(i==0);
    return extent_;
  }

  native_expr_t const & native() const{
    return nativeExprObj_;
  }

  native_expr_t & native(){
    return nativeExprObj_;
  }

  template<typename _MatrixType = MatrixType>
  mpl::enable_if_t<
    std::is_same<typename mytraits::memory_space, Kokkos::HostSpace>::value, 
    ref_t
    >
  operator()(size_t i) const
  {
    assert(i < (size_t)extent_);
    return nativeExprObj_(i);
  }
};
#endif


// #ifdef PRESSIO_ENABLE_TPL_PYBIND11
// template <typename MatrixType>
// struct DiagExpr<
//   MatrixType,
//   ::pressio::mpl::enable_if_t<
//     ::pressio::is_rank2_tensor_pybind<MatrixType>::value
//     >
//   >
// {
//   using this_t = DiagExpr<MatrixType>;
//   using traits = typename traits<this_t>;
//   using sc_t = typename traits::scalar_t;
//   using size_t = typename traits::size_t;
//   using ref_t = typename traits::reference_t;
//   using const_ref_t = typename traits::const_reference_t;
//   using pair_t = std::pair<std::size_t, std::size_t>;

// private:
//   std::reference_wrapper<MatrixType> matObj_;
//   size_t extent_ = {};

// public:
//   DiagExpr() = delete;
//   DiagExpr(const DiagExpr & other) = default;
//   DiagExpr & operator=(const DiagExpr & other) = delete;
//   DiagExpr(DiagExpr && other) = default;
//   DiagExpr & operator=(DiagExpr && other) = delete;
//   ~DiagExpr() = default;

//   DiagExpr(MatrixType & matObjIn)
//     : matObj_(matObjIn),
//       extent_(matObjIn.extent(0))
//   {
//     assert(matObjIn.extent(0) == matObjIn.extent(1));
//   }

// public:
//   size_t extent() const{
//     return extent_;
//   }

//   size_t extent(size_t i) const{
//     assert(i==0);
//     return extent_;
//   }

//   // non-const subscripting
//   template<typename _MatrixType = MatrixType>
//   mpl::enable_if_t<
//     !std::is_const<typename std::remove_reference<_MatrixType>::type>::value,
//     ref_t
//     >
//   operator()(size_t i)
//   {
//     assert(i < (size_t)extent_);
//     return matObj_.get()(i,i);
//   }

//   // const subscripting
//   const_ref_t operator()(size_t i) const
//   {
//     assert(i < (size_t)extent_);
//     return matObj_.get()(i,i);
//   }
// };
// #endif

}}}
#endif  // CONTAINERS_EXPRESSIONS_DIAG_CONTAINERS_DIAG_CLASSES_HPP_
