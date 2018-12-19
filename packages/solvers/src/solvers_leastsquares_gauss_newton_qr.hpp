
#ifndef SOLVERS_LEASTSQUARES_GAUSS_NEWTON_QR_HPP
#define SOLVERS_LEASTSQUARES_GAUSS_NEWTON_QR_HPP

#include "solvers_ConfigDefs.hpp"
#include "./experimental/solvers_meta_static_checks.hpp"
#include "../../CORE_OPS"
#include "../../core/src/meta/core_meta_detection_idiom.hpp"
#include "../../QR_BASIC"

namespace rompp{ namespace solvers{


template <typename scalar_t,
	  typename qr_type,
	  core::meta::enable_if_t<
	    core::meta::is_default_constructible<qr_type>::value
	    > * = nullptr
	  >
class GaussNewtonQR {

  using uint_t = core::default_types::uint;
  // using eig_col_vec = Eigen::Matrix<scalar_t, QR_n, 1>;
  // using eig_R_mat = Eigen::Matrix<scalar_t, QR_n, QR_n>;
  // using R_type = rompp::core::Matrix<eig_R_mat>;
  // template <typename jac_t>
  // using qr_t = rompp::qr::QRSolver<jac_t,
  // 				   rompp::core::MultiVector,
  // 				   R_type,
  // 				   qr_algo_tag>;

private:
  //core::Vector<eig_col_vec> QTResid_ = {};
  uint_t maxNonLinearIterations_     = 500;
  scalar_t nonLinearTolerance_       = 1e-6;
  scalar_t normO_ 		     = {};
  scalar_t normN_ 		     = {};
  qr_type qrObj 		     = {};

public:
  GaussNewtonQR() = default;
  GaussNewtonQR(const GaussNewtonQR &) = delete;
  ~GaussNewtonQR() = default;

public:
  void setMaxNonLinearIterations(uint_t maxNonLinearIterations) {
    maxNonLinearIterations_ = maxNonLinearIterations;
  }

  void setNonLinearTolerance(scalar_t nonLinearTolerance) {
    nonLinearTolerance_ = std::abs(nonLinearTolerance);
  }

  template <typename SystemT,
	    typename VectorT,
	    core::meta::enable_if_t<
	      core::meta::is_core_vector_wrapper<VectorT>::value and
	      std::is_same<VectorT, typename SystemT::vector_type >::value
	      > * =nullptr
	    >
  void solve(const SystemT& sys, VectorT& x)
  {
#ifdef DEBUG_PRINT
    std::cout << " starting Gauss-Newton solve "
	      << " tol = " << nonLinearTolerance_
	      << " maxIter = " << maxNonLinearIterations_
	      << std::endl;
#endif

    auto Residual = sys.residual(x);
    auto Jac = sys.jacobian(x);
    //using jac_t = decltype(Jac);
    // // create qr object
    // qr_t<jac_t> qrObj;

    VectorT QTResid_;
    auto dx(x);
    normO_ = ::rompp::core::ops::norm2(dx);
    normN_ = static_cast<scalar_t>(0);

    uint_t iStep = 1;
    while (iStep++ < maxNonLinearIterations_)
    {
      // // QR decomposition of Jacobian
      // qrObj.computeThin(Jac);
      // qrObj.project(QTResid_);
      // qrObj.solve(*dx.data(), QTResid_);

      // const auto & QF = qrObj.cRefQFactor();
      // const auto & RF = qrObj.cRefRFactor();
      // // extract Rn block from R factor
      // auto n = RF.cols();
      // auto RFn = RF.data()->block(0,0,n,n);

      // //----------------------------
      // // Solve: RFn dx = (Q^T Res)n
      // //----------------------------
      // // compute: Q^T * Residual
      // core::ops::dot(QF, Residual, QTResid_);
      // // get the n block of RHS
      // auto RHSn = QTResid_.data()->block(0,0,n,1);
      // *dx.data() = RFn.template triangularView<Eigen::Upper>().solve(RHSn);

      // update solution
      x -= dx;

      normN_ = ::rompp::core::ops::norm2(dx);
#ifdef DEBUG_PRINT
      std::cout << " GN step=" << iStep
		<< " norm(dx)= " << normN_
		<< std::endl;
#endif

      if (std::abs(normO_ - normN_) < nonLinearTolerance_){
#ifdef DEBUG_PRINT
	std::cout << " GN converged! " <<
		  << " final norm(dx)= " << normN_
		  << std::endl;
#endif
      	break;
      }
      normO_ = normN_;
      sys.residual(x, Residual);
      sys.jacobian(x, Jac);
    }

  }//solve

};//class

}} //end namespace rompp::solvers
#endif
