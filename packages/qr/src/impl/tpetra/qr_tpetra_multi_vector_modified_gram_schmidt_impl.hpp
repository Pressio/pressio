/*
//@HEADER
// ************************************************************************
//
// qr_tpetra_multi_vector_modified_gram_schmidt_impl.hpp
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

#ifndef QR_IMPL_TPETRA_QR_TPETRA_MULTI_VECTOR_MODIFIED_GRAM_SCHMIDT_IMPL_HPP_
#define QR_IMPL_TPETRA_QR_TPETRA_MULTI_VECTOR_MODIFIED_GRAM_SCHMIDT_IMPL_HPP_

namespace pressio{ namespace qr{ namespace impl{

template<typename matrix_t, typename R_t, typename MV_t, template<typename...> class Q_type>
class ModGramSchmidtMVTpetra<matrix_t, R_t, MV_t, Q_type, void>
{

  using this_t	     = ModGramSchmidtMVTpetra<matrix_t, R_t, MV_t, Q_type, void>;
  using int_t	     = int;
  using sc_t	     = typename containers::details::traits<matrix_t>::scalar_t;
  using eig_dyn_mat  = Eigen::Matrix<sc_t, Eigen::Dynamic, Eigen::Dynamic>;
  using R_nat_t	     = containers::DenseMatrix<eig_dyn_mat>;
  using Q_t	     = Q_type<MV_t>;

public:
  ModGramSchmidtMVTpetra() = default;
  ~ModGramSchmidtMVTpetra() = default;

  void computeThinOutOfPlace(const matrix_t & Ain)
  {
    auto & A = const_cast<matrix_t &>(Ain);

    auto nVecs = A.numVectors();
    auto & ArowMap = *A.data()->getMap();
    createQIfNeeded(ArowMap, nVecs);
    createLocalRIfNeeded(nVecs);

    sc_t rkkInv = {};
    for (auto k=0; k<A.numVectors(); k++)
    {
      auto ak = A.data()->getVector(k);
      localR_(k,k) = ak->norm2();
      rkkInv = utils::constants<sc_t>::one()/localR_(k,k);

      auto qk = Qmat_->data()->getVectorNonConst(k);
      qk->update( rkkInv, *ak, utils::constants<sc_t>::zero() );

      for (auto j=k+1; j<A.numVectors(); j++){
      	auto aj = A.data()->getVectorNonConst(j);
      	localR_(k,j) = qk->dot(*aj);
      	aj->update(-localR_(k,j), *qk, utils::constants<sc_t>::one());
      }
    }
  }

  template <typename vector_t>
  void doLinSolve(const vector_t & rhs, vector_t & y)const {
    //auto vecSize = y.size();
    auto & Rm = localR_.data()->template triangularView<Eigen::Upper>();
    *y.data() = Rm.solve(*rhs.data());
  }

  template < typename vector_in_t, typename vector_out_t>
  void applyQTranspose(const vector_in_t & vecIn, vector_out_t & vecOut) const
  {
    constexpr auto beta  = ::pressio::utils::constants<sc_t>::zero();
    constexpr auto alpha = ::pressio::utils::constants<sc_t>::one();
    ::pressio::ops::product(::pressio::transpose(), alpha, *this->Qmat_, vecIn, beta, vecOut);
  }

  const Q_t & QFactor() const {
    return *this->Qmat_;
  }

private:
  void createLocalRIfNeeded(int newsize){
    if (localR_.extent(0)!=newsize or localR_.extent(1)!=newsize){
      localR_ = R_nat_t(newsize, newsize);
      ::pressio::ops::set_zero(localR_);
    }
  }

  template <typename map_t>
  void createQIfNeeded(const map_t & map, int cols){
    if (!Qmat_ or !Qmat_->data()->getMap()->isSameAs(map) )
      Qmat_ = std::make_shared<Q_t>(map, cols);
  }

private:
  R_nat_t localR_			= {};
  // todo: these must be moved somewhere else
  mutable std::shared_ptr<Q_t> Qmat_	= nullptr;
  mutable std::shared_ptr<R_t> Rmat_	= nullptr;
};

}}} // end namespace pressio::qr::impl
#endif  // QR_IMPL_TPETRA_QR_TPETRA_MULTI_VECTOR_MODIFIED_GRAM_SCHMIDT_IMPL_HPP_
