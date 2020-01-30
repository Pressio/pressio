#include <gtest/gtest.h>
#include "CONTAINERS_ALL"
#include "ODE_ALL"
#include "SOLVERS_NONLINEAR"
#include "ROM_WLS"

struct ValidApp{
  using eigVec = Eigen::VectorXd;
  using scalar_type = double;
  using state_type  = eigVec;
  using velocity_type = eigVec;
  using dense_matrix_type = Eigen::MatrixXd;
  using jacobian_type = Eigen::MatrixXd;

public:
  void velocity(const state_type & y, scalar_type t, velocity_type & f) const
  {};

  velocity_type velocity(const state_type & y, scalar_type t) const{
    velocity_type f;
    return f;
  };

  void applyJacobian(const state_type & y,
		     const dense_matrix_type & B,
		     scalar_type t,
		     dense_matrix_type & A) const
  {}

  dense_matrix_type applyJacobian(const state_type & y,
				  const dense_matrix_type & B,
				  scalar_type t) const{
    dense_matrix_type A;
    return A;
  }
};

TEST(rom_wls_meta, hessAndGradSystem){
  using fom_t    = ValidApp;
  using scalar_t     = typename fom_t::scalar_type;
  using fom_native_state_t = typename fom_t::state_type;
  using fom_state_t        = ::pressio::containers::Vector<fom_native_state_t>;

  using eig_dyn_mat    = Eigen::Matrix<scalar_t, -1, -1>;
  using eig_dyn_vec    = Eigen::Matrix<scalar_t, -1, 1>;

  using wls_state_t    = pressio::containers::Vector<eig_dyn_vec>;
  using hessian_t          = pressio::containers::Matrix<eig_dyn_mat>;
  using decoder_jac_t    = pressio::containers::MultiVector<eig_dyn_mat>;
  using decoder_t    = pressio::rom::LinearDecoder<decoder_jac_t, wls_state_t, fom_state_t>;

  using bdf1_tag      = ::pressio::ode::implicitmethods::Euler;
  using bdf2_tag      = ::pressio::ode::implicitmethods::BDF2;
  using explicitEuler_tag      = ::pressio::ode::explicitmethods::Euler;

  using wls_system_t1 = pressio::rom::wls::SystemHessianAndGradientApi<fom_t,wls_state_t,decoder_t,bdf1_tag,hessian_t>;
  using wls_system_t2 = pressio::rom::wls::SystemHessianAndGradientApi<fom_t,wls_state_t,decoder_t,bdf2_tag,hessian_t>;
  using wls_system_t3 = pressio::rom::wls::SystemHessianAndGradientApi<fom_t,wls_state_t,decoder_t,explicitEuler_tag,hessian_t>;

  fom_t appObj;
  fom_state_t yFOM_IC(2);
  fom_state_t yFOM_Ref(2);
  decoder_jac_t Phi(2,1);
  decoder_t decoderObj(Phi);
  constexpr int numStepsInWindow = 5;
  wls_system_t1 wlsSystem_IE(appObj, yFOM_IC, yFOM_Ref, decoderObj, numStepsInWindow);
  wls_system_t2 wlsSystem_BDF2(appObj, yFOM_IC, yFOM_Ref, decoderObj, numStepsInWindow);
  wls_system_t3 wlsSystem_EE(appObj, yFOM_IC, yFOM_Ref, decoderObj, numStepsInWindow);


}

