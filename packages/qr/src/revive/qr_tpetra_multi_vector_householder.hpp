/*
//@HEADER
// ************************************************************************
//
// qr_tpetra_multi_vector_householder.hpp
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
#ifndef QR_TPETRA_MULTI_VECTOR_HOUSEHOLDER_HPP_
#define QR_TPETRA_MULTI_VECTOR_HOUSEHOLDER_HPP_

#include "qr_ConfigDefs.hpp"
#include "qr_fwd.hpp"
#include "qr_solver_base.hpp"

#include "../../CONTAINERS_ALL"
#include <Eigen/OrderingMethods>
#include<Eigen/SparseQR>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>


namespace pressio{ namespace qr{

// Tpetra multivector, householder
template<typename matrix_type,
	 typename R_type,
	 template <typename...> class Q_type>
class QRSolver<matrix_type,
	       ::pressio::qr::Householder,
	       R_type,
	       Q_type,
	       typename
	       std::enable_if<
		 containers::meta::is_multi_vector_wrapper_tpetra<matrix_type>::value and
		 containers::meta::is_matrix_wrapper<R_type>::value and
		 containers::details::traits<R_type>::is_shared_mem and
		 containers::details::traits<R_type>::is_dense
		 >::type
	       >
  : public QRSolverBase<QRSolver<matrix_type, ::pressio::qr::Householder, R_type, Q_type>,
			R_type,
			Q_type<typename containers::details::traits<matrix_type>::wrapped_t>,
			matrix_type>{

  using MV = typename containers::details::traits<matrix_type>::wrapped_t;
  using sc_t = typename containers::details::traits<matrix_type>::scalar_t;
  using LO_t = typename containers::details::traits<matrix_type>::local_ordinal_t;
  using GO_t = typename containers::details::traits<matrix_type>::global_ordinal_t;
  using map_t = typename containers::details::traits<matrix_type>::data_map_t;
  using node_t = typename containers::details::traits<matrix_type>::node_t;
  using hexsp = typename containers::details::traits<matrix_type>::host_exec_space_t;
  using Q_t = Q_type<MV>;

  using this_t = QRSolver<matrix_type, ::pressio::qr::Householder, R_type, Q_type>;
  using base_t = QRSolverBase<this_t, R_type, Q_t, matrix_type>;

  using base_t::Qmat_;
  using base_t::Rmat_;

public:
  QRSolver() = default;
  ~QRSolver() = default;

private:

  void computeThinImpl(matrix_type & A){

    auto m = A.globalLength();
    auto n = A.globalNumVectors();
    auto ArowMap = A.getRCPDataMap();
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
      Teuchos::rcp (new Teuchos::MpiComm<int> (MPI_COMM_SELF));

    // convert it to replicated eptra matrix
    using local_map_t = Tpetra::Map<LO_t, GO_t, node_t>;
    using rcp_local_map_t = Teuchos::RCP<const local_map_t>;
    rcp_local_map_t rcp_local_map = Teuchos::rcp( new local_map_t(m, 0, comm) );

    using import_t = Tpetra::Import<LO_t, GO_t, node_t>;
    import_t importer(ArowMap, rcp_local_map);
    matrix_type A2(rcp_local_map, n);
    A2.data()->doImport(*A.data(), importer, Tpetra::INSERT);

    // store it into an Eigen matrix
    containers::Matrix<Eigen::MatrixXd> eA2W(m,n);
    for (int j=0;j<n;j++){
      auto colData = A2.data()->getData(j);
      for (int i=0;i<m;i++)
    	eA2W(i,j) = colData[i];
    }

    // do QR in Eigen
    Eigen::HouseholderQR<Eigen::MatrixXd> eQR(*eA2W.data());
    auto Qm = eQR.householderQ() * Eigen::MatrixXd::Identity(m,n);
    auto & Rm = eQR.matrixQR().template triangularView<Eigen::Upper>();

    // store R factor
    //auto RFn = Rm.block(0,0,n,n);
    Rmat_ = std::make_shared<R_type>(Rm);

    // store Q into replicated Tpetra::Multivector
    Q_t locQ( rcp_local_map, Qm.cols() );
    auto trilD = locQ.data();
    trilD->template sync<Kokkos::HostSpace>();

    auto v2d = trilD->template getLocalView<Kokkos::HostSpace>();
    auto c0 = Kokkos::subview(v2d, Kokkos::ALL(), 0);
    // //we are going to change the host view
    trilD->template modify<Kokkos::HostSpace>();
    for (int i=0;i<Qm.rows();i++)
      for (int j=0;j<Qm.cols();j++)
    	v2d(i,j) = Qm(i,j);

    // import from local to distributed
    Qmat_ = std::make_shared<Q_t>(ArowMap, Qm.cols());
    import_t importer2(rcp_local_map, ArowMap);
    Qmat_->data()->doImport(*locQ.data(), importer2, Tpetra::INSERT);

  }

  const Q_t & cRefQFactorImpl() const {
    return *Qmat_;
  }

  const R_type & cRefRFactorImpl() const {
    return *Rmat_;
  }

  friend base_t;

};//end class

}} // end namespace pressio::qr
#endif
#endif
