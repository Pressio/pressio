
#ifndef ODE_MODELS_FOR_STATIC_CHECKS_FOR_IMPLICIT_STEPPER_HPP_
#define ODE_MODELS_FOR_STATIC_CHECKS_FOR_IMPLICIT_STEPPER_HPP_

namespace pressio{ namespace ode{ namespace testing{

struct ModelForImplicitMissingScalarTypedef{
  using scalar_t = double;
  using state_type    = Eigen::VectorXd;
  using velocity_type = Eigen::VectorXd;
  using jacobian_type = Eigen::MatrixXd;
public:
  void velocity(const state_type &, scalar_t , velocity_type &) const{}
  velocity_type velocity(const state_type &, scalar_t) const{ return velocity_type(); }
  void jacobian(const state_type &, scalar_t, jacobian_type &) const{}
  jacobian_type jacobian(const state_type &, scalar_t) const{ return jacobian_type(); }
};

struct ModelForImplicitMissingStateTypedef{
  using scalar_type = double;
  using state_t    = Eigen::VectorXd;
  using velocity_type = Eigen::VectorXd;
  using jacobian_type = Eigen::MatrixXd;
public:
  void velocity(const state_t &, scalar_type , velocity_type &) const{}
  velocity_type velocity(const state_t &, scalar_type) const{ return velocity_type(); }
  void jacobian(const state_t &, scalar_type, jacobian_type &) const{}
  jacobian_type jacobian(const state_t &, scalar_type) const{ return jacobian_type(); }
};

struct ModelForImplicitMissingVelocityTypedef{
  using scalar_type = double;
  using state_type    = Eigen::VectorXd;
  using velocity_t = Eigen::VectorXd;
  using jacobian_type = Eigen::MatrixXd;
public:
  void velocity(const state_type &, scalar_type , velocity_t &) const{}
  velocity_t velocity(const state_type &, scalar_type) const{ return velocity_t(); }
  void jacobian(const state_type &, scalar_type, jacobian_type &) const{}
  jacobian_type jacobian(const state_type &, scalar_type) const{ return jacobian_type(); }
};

struct ModelForImplicitMissingJacobianTypedef{
  using scalar_type = double;
  using state_type    = Eigen::VectorXd;
  using velocity_type = Eigen::VectorXd;
  using jacobian_t = Eigen::MatrixXd;
public:
  void velocity(const state_type &, scalar_type , velocity_type &) const{}
  velocity_type velocity(const state_type &, scalar_type) const{ return velocity_type(); }
  void jacobian(const state_type &, scalar_type, jacobian_t &) const{}
  jacobian_t jacobian(const state_type &, scalar_type) const{ return jacobian_t(); }
};

struct ModelForImplicitMissingVelocityMethods{
  using scalar_type = double;
  using state_type    = Eigen::VectorXd;
  using velocity_type = Eigen::VectorXd;
  using jacobian_type = Eigen::MatrixXd;
public:
  void jacobian(const state_type &, scalar_type, jacobian_type &) const{}
  jacobian_type jacobian(const state_type &, scalar_type) const{ return jacobian_type(); }
};

struct ModelForImplicitMissingVelocityNonVoidMethod{
  using scalar_type = double;
  using state_type    = Eigen::VectorXd;
  using velocity_type = Eigen::VectorXd;
  using jacobian_type = Eigen::MatrixXd;
public:
  void velocity(const state_type &, scalar_type , velocity_type &) const{}
  void jacobian(const state_type &, scalar_type, jacobian_type &) const{}
  jacobian_type jacobian(const state_type &, scalar_type) const{ return jacobian_type(); }
};

struct ModelForImplicitMissingVelocityVoidMethod{
  using scalar_type = double;
  using state_type    = Eigen::VectorXd;
  using velocity_type = Eigen::VectorXd;
  using jacobian_type = Eigen::MatrixXd;
public:
  velocity_type velocity(const state_type &, scalar_type) const{ return velocity_type(); }
  void jacobian(const state_type &, scalar_type, jacobian_type &) const{}
  jacobian_type jacobian(const state_type &, scalar_type) const{ return jacobian_type(); }
};

struct ModelForImplicitMissingJacobianMethods{
  using scalar_type = double;
  using state_type    = Eigen::VectorXd;
  using velocity_type = Eigen::VectorXd;
  using jacobian_type = Eigen::MatrixXd;
public:
  void velocity(const state_type &, scalar_type , velocity_type &) const{}
  velocity_type velocity(const state_type &, scalar_type) const{ return velocity_type(); }
};

struct ModelForImplicitMissingJacobianNonVoidMethod{
  using scalar_type = double;
  using state_type    = Eigen::VectorXd;
  using velocity_type = Eigen::VectorXd;
  using jacobian_type = Eigen::MatrixXd;
public:
  void velocity(const state_type &, scalar_type , velocity_type &) const{}
  velocity_type velocity(const state_type &, scalar_type) const{ return velocity_type(); }
  void jacobian(const state_type &, scalar_type, jacobian_type &) const{}
};

struct ModelForImplicitMissingJacobianVoidMethod{
  using scalar_type = double;
  using state_type    = Eigen::VectorXd;
  using velocity_type = Eigen::VectorXd;
  using jacobian_type = Eigen::MatrixXd;
public:
  void velocity(const state_type &, scalar_type , velocity_type &) const{}
  velocity_type velocity(const state_type &, scalar_type) const{ return velocity_type(); }
  jacobian_type jacobian(const state_type &, scalar_type) const{ return jacobian_type(); }
};

struct ModelForImplicitValid{
  using scalar_type   = double;
  using state_type    = Eigen::VectorXd;
  using velocity_type = Eigen::VectorXd;
  using jacobian_type = Eigen::MatrixXd;
public:
  void velocity(const state_type &, scalar_type , velocity_type &) const{}
  velocity_type velocity(const state_type &, scalar_type) const{ return velocity_type(); }
  void jacobian(const state_type &, scalar_type, jacobian_type &) const{}
  jacobian_type jacobian(const state_type &, scalar_type) const{ return jacobian_type(); }
};

}}} // namespace pressio::ode::testing
#endif
