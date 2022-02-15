/*
//@HEADER
// ************************************************************************
//
// qr_tpetra_mv_householder_using_eigen_impl.hpp
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

#ifndef QR_IMPL_TPETRA_QR_TPETRA_MV_HOUSEHOLDER_USING_EIGEN_IMPL_HPP_
#define QR_IMPL_TPETRA_QR_TPETRA_MV_HOUSEHOLDER_USING_EIGEN_IMPL_HPP_

#include <Eigen/OrderingMethods>
#include <Eigen/SparseQR>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>

namespace pressio{ namespace qr{ namespace impl{

template<typename MatrixType, typename R_t>
class TpetraMVHouseholderUsingEigen
{

public:
  using sc_t = typename ::pressio::Traits<MatrixType>::scalar_type;
  using lo_t         = typename ::pressio::Traits<MatrixType>::local_ordinal_type;
  using go_t         = typename ::pressio::Traits<MatrixType>::global_ordinal_type;
  using node_t       = typename ::pressio::Traits<MatrixType>::node_type;
  using Q_type       = Tpetra::MultiVector<sc_t, lo_t, go_t, node_t>;

  using eig_dyn_mat	= Eigen::Matrix<sc_t, -1, -1>;
  using help_impl_t	= QRHouseholderDenseEigenMatrix<eig_dyn_mat, R_t>;
  help_impl_t myImpl_	= {};

public:
  TpetraMVHouseholderUsingEigen() = default;
  ~TpetraMVHouseholderUsingEigen() = default;

  template < typename vector_in_t, typename vector_out_t>
  void applyQTranspose(const vector_in_t & vecIn, vector_out_t & vecOut) const
  {
    constexpr auto beta  = ::pressio::utils::Constants<sc_t>::zero();
    constexpr auto alpha = ::pressio::utils::Constants<sc_t>::one();
    ::pressio::ops::product(::pressio::transpose(), alpha, *this->Qmat_, vecIn, beta, vecOut);
  }

  template <typename VectorType>
  void doLinSolve(const VectorType & rhs, VectorType & y)const{
    myImpl_.template doLinSolve<VectorType>(rhs, y);
  }

  const Q_type & QFactor() const {
    return *this->Qmat_;
  }

  void computeThinOutOfPlace(const MatrixType & Ain)
  {
    auto & A = const_cast<MatrixType &>(Ain);

    auto rows = ::pressio::ops::extent(A,0);
    auto cols = ::pressio::ops::extent(A,1);
    auto ArowMap = A.getMap();
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
      Teuchos::rcp (new Teuchos::MpiComm<int> (MPI_COMM_SELF));

    // convert it to replicated eptra matrix
    using local_map_t = Tpetra::Map<lo_t, go_t, node_t>;
    using rcp_local_map_t = Teuchos::RCP<const local_map_t>;
    rcp_local_map_t rcp_local_map = Teuchos::rcp( new local_map_t(rows, 0, comm) );

    using import_t = Tpetra::Import<lo_t, go_t, node_t>;
    import_t importer(ArowMap, rcp_local_map);
    MatrixType A2(rcp_local_map, cols);
    A2.doImport(A, importer, Tpetra::INSERT);

    // store it into an Eigen matrix
    Eigen::Matrix<sc_t, -1, -1> eA2W(rows,cols);
    for (int j=0;j<cols;j++)
    {
      auto colData = A2.getData(j);
      for (int i=0;i<rows;i++){
    	 eA2W(i,j) = colData[i];
      }
    }

    myImpl_.computeThinOutOfPlace(eA2W);

    // store Q into replicated Tpetra::Multivector
    const auto & Q2 = myImpl_.QFactor();
    Q_type locQ( rcp_local_map, Q2.cols() );
    // // auto trilD = locQ.data();
    // locQ.template sync<Kokkos::HostSpace>();

    auto v2d = locQ.template getLocalView<Kokkos::HostSpace>(Tpetra::Access::ReadWriteStruct());
    auto c0 = Kokkos::subview(v2d, Kokkos::ALL(), 0);
    // // //we are going to change the host view
    // locQ.template modify<Kokkos::HostSpace>();
    for (int i=0;i<Q2.rows();i++)
      for (int j=0;j<Q2.cols();j++)
    	v2d(i,j) = Q2(i,j);

    // import from local to distributed
    Qmat_ = std::make_shared<Q_type>(ArowMap, Q2.cols());
    import_t importer2(rcp_local_map, ArowMap);
    Qmat_->doImport(locQ, importer2, Tpetra::INSERT);
  }

private:
  // todo: these must be moved somewhere else
  mutable std::shared_ptr<Q_type> Qmat_	= nullptr;
  mutable std::shared_ptr<R_t> Rmat_	= nullptr;

};//end class

}}} // end namespace pressio::qr::impl
#endif  // QR_IMPL_TPETRA_QR_TPETRA_MV_HOUSEHOLDER_USING_EIGEN_IMPL_HPP_
