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
#include <Epetra_Map.h>
#include <Epetra_LocalMap.h>

namespace pressio{ namespace qr{ namespace impl{

template<typename MatrixType, typename R_t>
class EpetraMVHouseholderUsingEigen
{

public:
  using sc_t = typename ::pressio::Traits<MatrixType>::scalar_type;
  using Q_type = Epetra_MultiVector;
  using eig_dyn_mat	= Eigen::MatrixXd;
  using help_impl_t	= QRHouseholderDenseEigenMatrix<eig_dyn_mat, R_t>;

private:
  help_impl_t myImpl_	= {};

public:
  EpetraMVHouseholderUsingEigen() = default;
  ~EpetraMVHouseholderUsingEigen() = default;

  template < typename VectorInType, typename VectorOutType>
  void applyQTranspose(const VectorInType & vecIn, VectorOutType & vecOut) const
  {
    constexpr auto beta  = ::pressio::utils::Constants<sc_t>::zero();
    constexpr auto alpha = ::pressio::utils::Constants<sc_t>::one();
    ::pressio::ops::product(::pressio::transpose(), alpha, *this->Qmat_, vecIn, beta, vecOut);
  }

  template <typename VectorType>
  void doLinSolve(const VectorType & rhs, VectorType & y)const
  {
    myImpl_.template doLinSolve<VectorType>(rhs, y);
  }

  const Q_type & QFactor() const 
  {
    return *this->Qmat_;
  }

  template < typename VectorInType, typename VectorOutType>
  void applyRTranspose(const VectorInType & vecIn, VectorOutType & y) const
  {
    myImpl_.applyRTranspose(vecIn, y);
  }

  void computeThinOutOfPlace(const MatrixType & A)
  {
    auto rows = ::pressio::ops::extent(A,0);
    auto cols = ::pressio::ops::extent(A,1);
    auto & ArowMap = A.Map();

    // convert it to replicated eptra matrix
    Epetra_LocalMap locMap(rows, 0, A.Comm());
    Epetra_Import importer(locMap, ArowMap);
    MatrixType A2(locMap, cols);
    A2.Import(A, importer, Insert);

    // store it into an Eigen matrix
    eig_dyn_mat eA2W(rows,cols);
    for (int i=0;i<rows;i++){
      for (int j=0;j<cols;j++){
    	 eA2W(i,j) = A2[j][i];
      }
    }

    myImpl_.computeThinOutOfPlace(eA2W);

    // store Q into replicated Epetra_Multivector
    const auto & Q2 = myImpl_.QFactor();
    Q_type locQ(locMap,Q2.cols());
    for (int i=0;i<Q2.rows();i++){
      for (int j=0;j<Q2.cols();j++){
    	 locQ[j][i] = Q2(i,j);
      }
    }

    // import from local to distributed
    Qmat_ = std::make_shared<Q_type>(ArowMap, Q2.cols());
    Epetra_Import importer2(ArowMap, locMap);
    Qmat_->Import(locQ, importer2, Insert);
  }

private:
  mutable std::shared_ptr<Q_type> Qmat_	= nullptr;
  mutable std::shared_ptr<R_t> Rmat_	= nullptr;
};

}}} // end namespace pressio::qr::impl
#endif  // QR_IMPL_EPETRA_QR_EPETRA_MV_HOUSEHOLDER_USING_EIGEN_IMPL_HPP_
