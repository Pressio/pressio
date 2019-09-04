/*
//@HEADER
// ************************************************************************
//
// qr_eigen_dense_matrix.hpp
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

#ifndef QR_EIGEN_DENSE_MATRIX_HPP_
#define QR_EIGEN_DENSE_MATRIX_HPP_

#include "qr_ConfigDefs.hpp"
#include "qr_fwd.hpp"
#include "qr_solver_base.hpp"
// #include "../../CONTAINERS_ALL"
#include <Eigen/OrderingMethods>
#include<Eigen/SparseQR>

namespace pressio{ namespace qr{

// the input data is a wrapper for eigen DENSE matrix
template<typename matrix_type,
	 typename R_type,
	 template <typename...> class Q_type>
class QRSolver<matrix_type,
	       ::pressio::qr::Householder,
	       R_type,
	       Q_type,
	       ::pressio::mpl::enable_if_t<
		 containers::meta::is_dense_matrix_wrapper_eigen<matrix_type>::value and
		 containers::meta::is_matrix_wrapper<R_type>::value and
		 containers::details::traits<R_type>::is_shared_mem and
		 containers::details::traits<R_type>::is_dense
		 >
	       >
  : public QRSolverBase<QRSolver<matrix_type, ::pressio::qr::Householder, R_type, Q_type>,
			R_type, Q_type<Eigen::MatrixXd>, matrix_type>{

  using sc_t = typename containers::details::traits<matrix_type>::scalar_t;
  using Q_t = Q_type<Eigen::MatrixXd>;
  using this_t = QRSolver<matrix_type, ::pressio::qr::Householder, R_type, Q_type>;
  using base_t = QRSolverBase<this_t, R_type, Q_t, matrix_type>;
  friend base_t;
  using base_t::Qmat_;
  using base_t::Rmat_;

public:
  QRSolver() = default;
  ~QRSolver() = default;

private:

  void computeThinImpl(matrix_type & A){

    auto m = A.rows();
    auto n = A.cols();

    using native_mat_type = typename containers::details::traits<matrix_type>::wrapped_t;
    Eigen::HouseholderQR<native_mat_type> eQR(*A.data());

    auto Qm = eQR.householderQ() * Eigen::MatrixXd::Identity(m,n);
    Qmat_ = std::make_shared<Q_t>(Qm);

    auto & Rm = eQR.matrixQR().template triangularView<Eigen::Upper>();
    Rmat_ = std::make_shared<R_type>( Rm );
  }

  const Q_t & cRefQFactorImpl() const {
    return *Qmat_;
  }

  const R_type & cRefRFactorImpl() const {
    return *Rmat_;
  }

};//end class

}} // end namespace pressio::qr
#endif
