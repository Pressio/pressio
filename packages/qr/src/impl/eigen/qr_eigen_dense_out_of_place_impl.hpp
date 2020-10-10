/*
//@HEADER
// ************************************************************************
//
// qr_eigen_dense_out_of_place_impl.hpp
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

#ifndef QR_IMPL_EIGEN_QR_EIGEN_DENSE_OUT_OF_PLACE_IMPL_HPP_
#define QR_IMPL_EIGEN_QR_EIGEN_DENSE_OUT_OF_PLACE_IMPL_HPP_

#include <Eigen/QR>

namespace pressio{ namespace qr{ namespace impl{

template< typename matrix_type, typename R_t, template <typename...> class Q_type>
class QRHouseholderDenseEigenMatrixWrapper<
  matrix_type, R_t, Q_type,
  ::pressio::mpl::enable_if_t<
    containers::predicates::is_dense_matrix_wrapper_eigen<matrix_type>::value or
    containers::predicates::is_multi_vector_wrapper_eigen<matrix_type>::value
    >
  >{

  using sc_t	     = typename containers::details::traits<matrix_type>::scalar_t;
  using nat_mat_t    = typename containers::details::traits<matrix_type>::wrapped_t;
  using factorizer_t = Eigen::HouseholderQR<nat_mat_t>;
  using Q_nat_t	     = Eigen::Matrix<sc_t, Eigen::Dynamic, Eigen::Dynamic>;
  using Q_t	     = Q_type<Q_nat_t>;

  mutable std::shared_ptr<Q_t> Qmat_	     = {};
  mutable std::shared_ptr<factorizer_t> fct_ = {};

public:
  QRHouseholderDenseEigenMatrixWrapper() = default;
  ~QRHouseholderDenseEigenMatrixWrapper() = default;

  void computeThinOutOfPlace(const matrix_type & A){
    auto rows = A.data()->rows();
    auto cols = A.data()->cols();
    fct_ = std::make_shared<factorizer_t>(*A.data());

    if (!Qmat_ or (Qmat_->data()->rows()!=rows and Qmat_->data()->cols()!=cols ) )
      Qmat_ = std::make_shared<Q_t>(rows,cols);

    *Qmat_->data() = fct_->householderQ() * Q_nat_t::Identity(rows,cols);
  }

  template < typename vector_in_t, typename vector_out_t>
  void applyQTranspose(const vector_in_t & vecIn, vector_out_t & vecOut) const
  {
    constexpr auto beta  = ::pressio::utils::constants<sc_t>::zero();
    constexpr auto alpha = ::pressio::utils::constants<sc_t>::one();
    ::pressio::ops::product(::pressio::transpose(), alpha, *this->Qmat_, vecIn, beta, vecOut);
  }

  template < typename vector_in_t, typename vector_out_t>
  void applyRTranspose(const vector_in_t & vecIn, vector_out_t & y) const
  {
    // y = R^T vecIn
    auto vecSize = y.extent(0);
    auto & Rm = fct_->matrixQR().block(0,0,vecSize,vecSize).template triangularView<Eigen::Upper>();
    *y.data() = Rm.transpose() * (*vecIn.data());
  }

  template <typename vector_t>
  void doLinSolve(const vector_t & rhs, vector_t & y)const{
    auto vecSize = y.extent(0);
    auto & Rm = fct_->matrixQR().block(0,0,vecSize,vecSize).
      template triangularView<Eigen::Upper>();
    *y.data() = Rm.solve(*rhs.data());
  }

  const Q_t & QFactor() const {
    return *this->Qmat_;
  }
};

}}} // end namespace pressio::qr::impl
#endif  // QR_IMPL_EIGEN_QR_EIGEN_DENSE_OUT_OF_PLACE_IMPL_HPP_
