/*
//@HEADER
// ************************************************************************
//
// qr_epetra_multi_vector_modified_gram_schmidt_impl.hpp
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

#ifndef QR_IMPL_EPETRA_QR_EPETRA_MULTI_VECTOR_MODIFIED_GRAM_SCHMIDT_IMPL_HPP_
#define QR_IMPL_EPETRA_QR_EPETRA_MULTI_VECTOR_MODIFIED_GRAM_SCHMIDT_IMPL_HPP_

namespace pressio{ namespace qr{ namespace impl{

template<typename MatrixType, typename R_t>
class ModGramSchmidtMVEpetra
{

public:
  using int_t	     = int;
  using sc_t	     = typename ::pressio::Traits<MatrixType>::scalar_type;
  using R_nat_t	    = Eigen::Matrix<sc_t, Eigen::Dynamic, Eigen::Dynamic>;
  using Q_type      = Epetra_MultiVector;
  static constexpr sc_t one_ = static_cast<sc_t>(1);
  static constexpr sc_t zero_ = static_cast<sc_t>(0);

public:
  ModGramSchmidtMVEpetra() = default;
  ~ModGramSchmidtMVEpetra() = default;

  void computeThinOutOfPlace(const MatrixType & Ain) 
  {
    auto & A = const_cast<MatrixType &>(Ain);

    auto nVecs = ::pressio::ops::extent(A,1);
    auto & ArowMap = A.Map();
    createQIfNeeded(ArowMap, nVecs);
    createLocalRIfNeeded(nVecs);

    sc_t rkkInv = zero_;
    for (auto k=0; k<nVecs; k++)
    {
      auto & ak = A(k);
      ak->Norm2(&localR_(k,k));
      rkkInv = one_/localR_(k,k);

      auto & qk = (*Qmat_)(k);
      qk->Update( rkkInv, *ak, zero_ );

      for (auto j=k+1; j<nVecs; j++){
	     auto & aj = A(j);
	     qk->Dot(*aj, &localR_(k,j));
	     aj->Update(-localR_(k,j), *qk, one_);
      }
    }
  }

  template <typename VectorType>
  void doLinSolve(const VectorType & rhs, VectorType & y)const {
    auto & Rm = localR_.template triangularView<Eigen::Upper>();
    y = Rm.solve(rhs);
  }

  template < typename VectorInType, typename VectorOutType>
  void applyQTranspose(const VectorInType & vecIn, VectorOutType & vecOut) const
  {
    constexpr auto beta  = ::pressio::utils::Constants<sc_t>::zero();
    constexpr auto alpha = ::pressio::utils::Constants<sc_t>::one();
    ::pressio::ops::product(::pressio::transpose(), alpha, *this->Qmat_, vecIn, beta, vecOut);
  }

  const Q_type & QFactor() const {
    return *this->Qmat_;
  }

private:
  void createLocalRIfNeeded(int newsize){
    const auto locRext0 = ::pressio::ops::extent(localR_, 0);
    const auto locRext1 = ::pressio::ops::extent(localR_, 1);   
    if (locRext0!=newsize or locRext1!=newsize){
      localR_ = R_nat_t(newsize, newsize);
      ::pressio::ops::set_zero(localR_);
    }
  }

  template <typename map_t>
  void createQIfNeeded(const map_t & map, int cols){
    if (!Qmat_ or !Qmat_->Map().SameAs(map))
      Qmat_ = std::make_shared<Q_type>(map, cols);
  }

private:
  R_nat_t localR_			= {};

  // todo: these must be moved somewhere else
  mutable std::shared_ptr<Q_type> Qmat_	= nullptr;
  mutable std::shared_ptr<R_t> Rmat_	= nullptr;
};

}}} // end namespace pressio::qr::impl
#endif  // QR_IMPL_EPETRA_QR_EPETRA_MULTI_VECTOR_MODIFIED_GRAM_SCHMIDT_IMPL_HPP_
