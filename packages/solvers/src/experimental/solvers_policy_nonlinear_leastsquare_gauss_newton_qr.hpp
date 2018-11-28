
#ifndef SOLVERS_EXPERIMENTAL_POLICY_NONLINEAR_LEASTSQUARE_GAUSS_NEWTON_QR_HPP
#define SOLVERS_EXPERIMENTAL_POLICY_NONLINEAR_LEASTSQUARE_GAUSS_NEWTON_QR_HPP

#include "solvers_system_traits.hpp"
#include "solvers_linear_factory.hpp"
#include "solvers_nonlinear_traits.hpp"
#include "solvers_l2_vector_norm.hpp"
#include "solvers_meta_static_checks.hpp"
#include "../solvers_ConfigDefs.hpp"
// #include "../../../CORE_BASIC"
#include "../../../CORE_OPS"
#include "../../../core/src/meta/core_meta_detection_idiom.hpp"
#include "../../../qr/src/qr_hacked.hpp"

namespace rompp{ namespace solvers{

    
struct SolversNonLinearIterativeLeastSquareGaussNewtonQRPolicy {

private:
  core::default_types::uint maxNonLinearIterations_ = 500;
  double nonLinearTolerance_ = 1e-6;


public:

  void setMaxNonLinearIterations(core::default_types::uint  maxNonLinearIterations) {
    maxNonLinearIterations_ = maxNonLinearIterations;
  }

  void setNonLinearTolerance(double nonLinearTolerance) {
    nonLinearTolerance_ = std::abs(nonLinearTolerance);
  }


  template <typename SystemT, 
	    typename VectorT,
	    typename core::meta::enable_if_t<
	      core::details::traits<VectorT>::is_vector &&
	      solvers::meta::are_vector_compatible<
		typename details::system_traits<SystemT>::vector_type,
		VectorT
		>::value
	      >* = nullptr
  >
  VectorT solve(const SystemT& sys,
		const VectorT& x0)
  {
    using sc_t = typename core::details::traits<VectorT>::scalar_t;
    using eig_mat = Eigen::Matrix<sc_t,Eigen::Dynamic,Eigen::Dynamic>;
    using eig_vec = Eigen::Matrix<sc_t,Eigen::Dynamic,1>;

    using NormT = L2Norm;
    
    auto Res = sys.residual(x0);
    auto Jac = sys.jacobian(x0);

    // for QR
    using jac_t = decltype(Jac);
    using R_type = rompp::core::Matrix<eig_mat>;
    rompp::qr::hack::QRSolver<jac_t, rompp::core::MultiVector, R_type> qrObj;
    core::Vector<eig_vec> QTRes;//(Jac.cols());
    
    auto x(x0);
    auto dx(x);
    double normN = 0.0;
    double normO = NormT::template compute_norm(dx);
    core::default_types::uint iStep = 1;
    while (iStep++ < maxNonLinearIterations_)
    { 
      // QR decomposition of Jacobian
      qrObj.compute(Jac);
      const auto & QF = qrObj.cRefQFactor();
      const auto & RF = qrObj.cRefRFactor();

      // extract Rn block from R factor
      auto n = RF.cols();
      auto RFn = RF.data()->block(0,0,n,n);
      
      // RFn dx = (Q^T Res)n
      core::ops::dot(QF, Res, QTRes); // compute: Q^T * Residual
      // get the n block of RHS
      auto RHSn = QTRes.data()->block(0,0,n,1);
      *dx.data() = RFn.template triangularView<Eigen::Upper>().solve(RHSn);

      // update solution 
      x -= dx;

      normN = NormT::template compute_norm(dx);
      // std::cout << " GN-iStep " << iStep 
      // 		<< " norm " << normN << " "
      // 		<< *x.data() << std::endl;

      if (std::abs(normO - normN) < nonLinearTolerance_){
      	break;
      }
      normO = normN;
      sys.residual(x, Res);
      sys.jacobian(x, Jac);
    }
    return x;
  }
};


}} // end namespace rompp::solvers
#endif
