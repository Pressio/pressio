/*
//@HEADER
// ************************************************************************
//
// qr_epetra_mv_householder_using_eigen_impl.hpp
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

#ifndef QR_IMPL_EPETRA_QR_EPETRA_MV_HOUSEHOLDER_USING_EIGEN_IMPL_HPP_
#define QR_IMPL_EPETRA_QR_EPETRA_MV_HOUSEHOLDER_USING_EIGEN_IMPL_HPP_

#include <Eigen/OrderingMethods>
#include <Eigen/SparseQR>
#include <Epetra_Import.h>

namespace pressio{ namespace qr{ namespace impl{

template<typename matrix_t, typename R_t, template <typename...> class Q_type>
class EpetraMVHouseholderUsingEigen
{

  using MV = Epetra_MultiVector;
  using sc_t = typename containers::details::traits<matrix_t>::scalar_t;
  using Q_t = Q_type<MV>;

  using eig_dyn_mat	= Eigen::MatrixXd;
  using eig_mat_w	= containers::DenseMatrix<eig_dyn_mat>;
  using help_impl_t	= QRHouseholderDenseEigenMatrixWrapper<eig_mat_w, R_t, Q_type>;
  help_impl_t myImpl_	= {};

public:
  EpetraMVHouseholderUsingEigen() = default;
  ~EpetraMVHouseholderUsingEigen() = default;

  template < typename vector_in_t, typename vector_out_t>
  void applyQTranspose(const vector_in_t & vecIn, vector_out_t & vecOut) const
  {
    constexpr auto beta  = ::pressio::utils::constants<sc_t>::zero();
    constexpr auto alpha = ::pressio::utils::constants<sc_t>::one();
    ::pressio::ops::product(::pressio::transpose(), alpha, *this->Qmat_, vecIn, beta, vecOut);
  }

  template <typename vector_t>
  void doLinSolve(const vector_t & rhs, vector_t & y)const{
    myImpl_.template doLinSolve<vector_t>(rhs, y);
  }

  const Q_t & QFactor() const {
    return *this->Qmat_;
  }

  template < typename vector_in_t, typename vector_out_t>
  void applyRTranspose(const vector_in_t & vecIn, vector_out_t & y) const
  {
    myImpl_.applyRTranspose(vecIn, y);
  }

  void computeThinOutOfPlace(const matrix_t & A)
  {
    auto rows = A.extent(0);
    auto cols = A.numVectors();
    auto & ArowMap = A.data()->Map();

    // convert it to replicated eptra matrix
    Epetra_LocalMap locMap(rows, 0, A.data()->Comm());
    Epetra_Import importer(locMap, ArowMap);
    matrix_t A2(locMap, cols);
    A2.data()->Import(*A.data(), importer, Insert);

    // store it into an Eigen matrix
    eig_mat_w eA2W(rows,cols);
    for (int i=0;i<rows;i++)
      for (int j=0;j<cols;j++)
    	eA2W(i,j) = A2(i,j);

    myImpl_.computeThinOutOfPlace(eA2W);

    // // do QR in Eigen
    // hhold_t eQR(*eA2W.data());
    // auto Qm = eQR.householderQ() * eig_dyn_mat::Identity(rows,n);
    // eig_dyn_mat Rm(eQR.matrixQR().template triangularView<Eigen::Upper>());

    // // store R factor
    // Rmat_ = std::make_shared<R_type>( Rm.block(0,0,n,n) );

    // store Q into replicated Epetra_Multivector
    const auto & Q2 = *(myImpl_.QFactor().data());
    Q_t locQ(locMap,Q2.cols());
    for (int i=0;i<Q2.rows();i++)
      for (int j=0;j<Q2.cols();j++)
    	locQ(i,j) = Q2(i,j);

    // import from local to distributed
    Qmat_ = std::make_shared<Q_t>(ArowMap, Q2.cols());
    Epetra_Import importer2(ArowMap, locMap);
    Qmat_->data()->Import(*locQ.data(), importer2, Insert);
  }

private:

  // todo: these must be moved somewhere else
  mutable std::shared_ptr<Q_t> Qmat_	= nullptr;
  mutable std::shared_ptr<R_t> Rmat_	= nullptr;

};//end class


}}} // end namespace pressio::qr::impl
#endif  // QR_IMPL_EPETRA_QR_EPETRA_MV_HOUSEHOLDER_USING_EIGEN_IMPL_HPP_
