/*
//@HEADER
// ************************************************************************
//
// qr_tpetra_block_multi_vector_tsqr_impl.hpp
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

#ifndef QR_IMPL_TPETRA_QR_TPETRA_BLOCK_MULTI_VECTOR_TSQR_IMPL_HPP_
#define QR_IMPL_TPETRA_QR_TPETRA_BLOCK_MULTI_VECTOR_TSQR_IMPL_HPP_

#include "Tpetra_TsqrAdaptor.hpp"

namespace pressio{ namespace qr{ namespace impl{

template<typename MatrixType, typename R_t>
class TpetraBlockMVTSQR
{

public:
  using int_t	     = int;
  using sc_t	     = typename ::pressio::Traits<MatrixType>::scalar_type;
  using serden_mat_t = Teuchos::SerialDenseMatrix<int_t, sc_t>;
  using trcp_mat     = Teuchos::RCP<serden_mat_t>;

  using lo_t	   = typename ::pressio::Traits<MatrixType>::local_ordinal_type;
  using go_t	   = typename ::pressio::Traits<MatrixType>::global_ordinal_type;
  using node_t   = typename ::pressio::Traits<MatrixType>::node_type;

  using Q_type	          = Tpetra::BlockMultiVector<sc_t, lo_t, go_t, node_t>;
  using tpetra_mv_t       = Tpetra::MultiVector<sc_t, lo_t, go_t, node_t>;
  using tsqr_adaptor_type = Tpetra::TsqrAdaptor<tpetra_mv_t>;

public:
  TpetraBlockMVTSQR() = default;
  ~TpetraBlockMVTSQR() = default;

  void computeThinOutOfPlace(const MatrixType & Ain)
  {
    auto & A = const_cast<MatrixType &>(Ain);

    auto nVecs	   = ::pressio::ops::extent(A,1);
    auto blockSize = A.getBlockSize();
    createLocalRIfNeeded(nVecs);

    // this is the row map of the block MV
    auto ArowMap = A.getMap();
    createQIfNeeded(*ArowMap, blockSize, nVecs);

    // get the multivector
    auto mv = A.getMultiVectorView();
    auto Qv = Qmat_->getMultiVectorView();
    tsqrAdaptor_.factorExplicit(mv, Qv, *localR_.get(), false);
  }

  template <typename VectorType>
  void doLinSolve(const VectorType & rhs, VectorType & y)const
  {
    qr::impl::solve<VectorType, trcp_mat>(rhs, this->localR_, y);
  }

  template < typename VectorInType, typename VectorOutType>
  void applyQTranspose(const VectorInType & vecIn, VectorOutType & vecOut) const
  {
    constexpr auto beta  = ::pressio::utils::Constants<sc_t>::zero();
    constexpr auto alpha = ::pressio::utils::Constants<sc_t>::one();
    ::pressio::ops::product(::pressio::transpose(), alpha, *this->Qmat_, vecIn, beta, vecOut);
  }

  template < typename VectorInType, typename VectorOutType>
  void applyRTranspose(const VectorInType & vecIn, VectorOutType & y) const
  {
    constexpr auto beta  = ::pressio::utils::Constants<sc_t>::zero();
    constexpr auto alpha = ::pressio::utils::Constants<sc_t>::one();
    ::pressio::ops::product(::pressio::transpose(), alpha, *this->localR_, vecIn, beta, y);
  }

  template <typename T = R_t>
  ::pressio::mpl::enable_if_t<
    !::pressio::is_dense_matrix_teuchos<T>::value and
    !std::is_void<T>::value,
    const T &
  >
  RFactor() const
  {
    this->Rmat_ = std::make_shared<T>(this->localR_->values());
    return *this->Rmat_;
  }

  template <typename T = R_t>
  ::pressio::mpl::enable_if_t<
    ::pressio::is_dense_matrix_teuchos<T>::value and
    !std::is_void<T>::value,
    const T &
  >
  RFactor() const
  {
    this->Rmat_ = std::make_shared<T>(*this->localR_, Teuchos::View);
    return *this->Rmat_;
  }

  const Q_type & QFactor() const {
    return *this->Qmat_;
  }

private:
  void createLocalRIfNeeded(int newsize)
  {
    if (localR_.is_null() or
    	(localR_->numRows()!=newsize and localR_->numCols()!=newsize)){
      localR_ = Teuchos::rcp(new serden_mat_t(newsize, newsize) );
    }
  }

  template <typename map_t>
  void createQIfNeeded(const map_t & map, int blockSize, int numVecs)
  {
    if (!Qmat_ or !Qmat_->getMap()->isSameAs(map) ){
     Qmat_ = std::make_shared<Q_type>(map, blockSize, numVecs);
    }
  }

private:
  tsqr_adaptor_type tsqrAdaptor_;
  trcp_mat localR_			= {};
  int computedRank_			= {};
  mutable std::shared_ptr<Q_type> Qmat_	= nullptr;
  mutable std::shared_ptr<R_t> Rmat_	= nullptr;
};

}}} // end namespace pressio::qr::impl
#endif  // QR_IMPL_TPETRA_QR_TPETRA_BLOCK_MULTI_VECTOR_TSQR_IMPL_HPP_
