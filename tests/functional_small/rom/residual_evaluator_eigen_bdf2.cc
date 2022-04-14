
#include <gtest/gtest.h>
#include "pressio/rom_decoder.hpp"
#include "pressio/rom_residual_evaluator.hpp"

struct TrivialFomContTimeEigen
{
  using scalar_type    = double;
  using state_type     = Eigen::VectorXd;
  using velocity_type  = state_type;
  int N_ = {};

  TrivialFomContTimeEigen(int N): N_(N){}

  velocity_type createVelocity() const{ return velocity_type(N_); }

  void velocity(const state_type & u, scalar_type time, velocity_type & f) const
  {
    EXPECT_TRUE((std::size_t)u.size()==(std::size_t)f.size());
    EXPECT_TRUE((std::size_t)u.size()==(std::size_t)N_);
    for (int i=0; i<f.rows(); ++i){
     f(i) = u(i) + time;
    }
  }
};

struct RomStateProducerBdf2
{
  using scalar_type = double;

  // M contains rom states for specific time steps e.g.
  // NOTE: what M contains is NOT necessarily
  // ALL time steps could be a subset

  // for BDF2 it could look like:
  //             n-2     n-1   n
  // M(:,0,2) = [y_t2,  y_t3,  y_t4]
  // M(:,3,5) = [y_t10, y_t11, y_t12]
  // M(:,6,8) = [y_t20, y_t21, y_t22]

  Eigen::MatrixXd M_;
  std::vector<double> times_;

  RomStateProducerBdf2(int romSize)
  {
    // times at n
    times_ = {10., 12., 15};
    M_.resize(romSize, 9);
    M_.col(0).fill(1.);
    M_.col(1).fill(2.);
    M_.col(2).fill(3.);
    M_.col(3).fill(4.);
    M_.col(4).fill(5.);
    M_.col(5).fill(6.);
    M_.col(6).fill(7.);
    M_.col(7).fill(8.);
    M_.col(8).fill(9.);
  }

  auto getMatrix() const { return M_; }

  auto getTimes()const { return times_; }

  // required
  std::size_t numberOfTimeSchemeStencilInstances() const {
    return 3;
  }

  // required
  scalar_type timeStepSizeAt(std::size_t stencilInstance) const
  {
         if (stencilInstance == 0){ return 0.1; }
    else if (stencilInstance == 1){ return 0.2; }
    else if (stencilInstance == 2){ return 0.3; }
    return {};
  }

  // required
  scalar_type timeAt(std::size_t stencilInstance, ::pressio::ode::n) const{
    return times_[stencilInstance];
  }

  scalar_type timeAt(std::size_t stencilInstance, ::pressio::ode::nMinusOne) const{
    return times_[stencilInstance] - timeStepSizeAt(stencilInstance);
  }

  // required
  auto stateAt(std::size_t stencilInstance, ::pressio::ode::n) const {
    return Eigen::VectorXd(M_.col(stencilInstance*3+2));
  }

  auto stateAt(std::size_t stencilInstance, ::pressio::ode::nMinusOne) const {
    return Eigen::VectorXd(M_.col(stencilInstance*3+1));
  }

  auto stateAt(std::size_t stencilInstance, ::pressio::ode::nMinusTwo) const {
    return Eigen::VectorXd(M_.col(stencilInstance*3));
  }
};

struct ObserverBdf2
{
  TrivialFomContTimeEigen fomObj_;
  Eigen::MatrixXd phi_;
  Eigen::MatrixXd M_;
  std::vector<double> times_;
  Eigen::VectorXd f_;
  Eigen::VectorXd R_gold_;

  ObserverBdf2(Eigen::MatrixXd phiIn,
	       Eigen::MatrixXd Min,
	       std::vector<double> times)
    : fomObj_(phiIn.rows()),
      phi_(phiIn),
      M_(Min),
      times_(times),
      f_(phiIn.rows()),
      R_gold_(phiIn.rows()){}

  template<class T>
  void operator()(std::size_t instanceId, const T & R)
  {
    if (instanceId==0){
      const auto ynm2 = phi_*M_.col(0);
      const auto ynm1 = phi_*M_.col(1);
      const auto yn   = phi_*M_.col(2);
      const auto timen = times_[0];
      fomObj_.velocity(yn, timen, f_);
      R_gold_ = yn - (4./3.)*ynm1 + (1./3.)*ynm2 - (2./3.)*0.1*f_;
      EXPECT_TRUE(R_gold_.isApprox(R));
    }

    if (instanceId==1){
      const auto ynm2 = phi_*M_.col(3);
      const auto ynm1 = phi_*M_.col(4);
      const auto yn   = phi_*M_.col(5);
      const auto timen = times_[1];
      fomObj_.velocity(yn, timen, f_);
      R_gold_ = yn - (4./3.)*ynm1 + (1./3.)*ynm2 - (2./3.)*0.2*f_;
      EXPECT_TRUE(R_gold_.isApprox(R));
    }

    if (instanceId==2){
      const auto ynm2 = phi_*M_.col(6);
      const auto ynm1 = phi_*M_.col(7);
      const auto yn   = phi_*M_.col(8);
      const auto timen = times_[2];
      fomObj_.velocity(yn, timen, f_);
      R_gold_ = yn - (4./3.)*ynm1 + (1./3.)*ynm2 - (2./3.)*0.3*f_;
      EXPECT_TRUE(R_gold_.isApprox(R));
    }
  }
};

TEST(rom, residual_evaluator_bdf2)
{
  using fom_state_t	  = Eigen::VectorXd;
  using decoder_jac_t	= Eigen::MatrixXd;

  int fomSize = 6;
  int romSize = 3;
  decoder_jac_t A(fomSize, romSize);
  A.col(0).fill(1.);
  A.col(1).fill(2.);
  A.col(2).fill(3.);

  namespace prom = pressio::rom;
  auto decoder = prom::create_time_invariant_linear_decoder<fom_state_t>(A);

  fom_state_t fomRefState(fomSize);
  fomRefState.fill(0.);

  auto scheme = ::pressio::ode::StepScheme::BDF2;
  auto eval = prom::create_residual_evaluator(scheme, decoder, fomRefState);

  TrivialFomContTimeEigen fomObj(fomSize);

  RomStateProducerBdf2 romP(romSize);
  auto R = fomObj.createVelocity();
  ObserverBdf2 obs(A, romP.getMatrix(), romP.getTimes());
  eval(fomObj, romP, R, obs);
}
