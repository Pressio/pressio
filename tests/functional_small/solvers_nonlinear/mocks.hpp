
#ifndef PRESSIO_FUNCTIONAL_SMALL_TESTS_NONLINEAR_SOLVERS_MOCKS_HPP_
#define PRESSIO_FUNCTIONAL_SMALL_TESTS_NONLINEAR_SOLVERS_MOCKS_HPP_

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "Eigen/SparseCore"
#include "pressio/solvers_nonlinear.hpp"

// =========================================
// problem class with residual-jacobian API 
// =========================================
template<bool use_dense_jacobian>
struct ConcreteSystemRJ
{
  using scalar_type = double;
  using state_type  = Eigen::VectorXd;
  using residual_type = state_type;
  using jacobian_type = typename std::conditional<
   use_dense_jacobian, 
   Eigen::MatrixXd, 
   Eigen::SparseMatrix<double>
   >::type;

 public:
  state_type createState() const { return residual_type(5); }
  residual_type createResidual() const { return residual_type(5); }
  jacobian_type createJacobian() const { return jacobian_type(5, 5); }
  void residual(const state_type& /*x*/, residual_type& res) const
  {
    res.setConstant(1);
  }

  template<bool _use_dense_jacobian = use_dense_jacobian>
  typename std::enable_if<!_use_dense_jacobian>::type 
  jacobian(const state_type& x, jacobian_type& jac) const
  {
    // (0  1 4  0  2)
    // (0  0 3 -2  2)
    // (1  1 0  0  2)
    // (4 -1 0  1  2)
    // (3  0 0  1  5)

    using trl_t = Eigen::Triplet<scalar_type>;
    std::vector<trl_t> tripletList;
    tripletList.push_back(trl_t(0,1,1.));
    tripletList.push_back(trl_t(0,2,4.));
    tripletList.push_back(trl_t(0,4,2.));
    tripletList.push_back(trl_t(1,2,3.));
    tripletList.push_back(trl_t(1,3,-2.));
    tripletList.push_back(trl_t(1,4,2.));
    tripletList.push_back(trl_t(2,0,1.));
    tripletList.push_back(trl_t(2,1,1.));
    tripletList.push_back(trl_t(2,4,2.));
    tripletList.push_back(trl_t(3,0,4.));
    tripletList.push_back(trl_t(3,1,-1.));
    tripletList.push_back(trl_t(3,3,1.));
    tripletList.push_back(trl_t(3,4,2.));
    tripletList.push_back(trl_t(4,0,3.));
    tripletList.push_back(trl_t(4,3,1.));
    tripletList.push_back(trl_t(4,4,5.));

    jac.setFromTriplets(tripletList.begin(), tripletList.end());
  }

  template<bool _use_dense_jacobian = use_dense_jacobian>
  typename std::enable_if<_use_dense_jacobian>::type 
  jacobian(const state_type& /*x*/, jacobian_type& jac) const
  {
    // (0  1 4  0  2)
    // (0  0 3 -2  2)
    // (1  1 0  0  2)
    // (4 -1 0  1  2)
    // (3  0 0  1  5)

    jac.setZero();
    jac(0,1) = 1.;
    jac(0,2) = 4.;
    jac(0,4) = 2.;
    jac(1,2) = 3.;
    jac(1,3) = -2.;
    jac(1,4) = 2.;
    jac(2,0) = 1.;
    jac(2,1) = 1.;
    jac(2,4) = 2.;
    jac(3,0) = 4.;
    jac(3,1) = -1.;
    jac(3,3) = 1.;
    jac(3,4) = 2.;
    jac(4,0) = 3.;
    jac(4,3) = 1.;
    jac(4,4) = 5.;
  }

  Eigen::MatrixXd goldHessian() const
  {
    auto J = createJacobian();
    state_type x;
    jacobian(x, J);
    Eigen::MatrixXd H = J.transpose() * J;
    return H;
  }
};

template<bool use_dense_jacobian>
struct MockSystemRJ
{
  using concrete_t = ConcreteSystemRJ<use_dense_jacobian>;
  using scalar_type   = typename concrete_t::scalar_type;
  using state_type    = typename concrete_t::state_type;
  using residual_type = typename concrete_t::residual_type;
  using jacobian_type = typename concrete_t::jacobian_type;

 public:
  MockSystemRJ()
  {
    ON_CALL(*this, createState).WillByDefault(
      [this]() {
      return mySystem_.createState();
    });
    ON_CALL(*this, createResidual).WillByDefault(
      [this]() {
      return mySystem_.createResidual();
    });
    ON_CALL(*this, createJacobian).WillByDefault(
      [this]() {
      return mySystem_.createJacobian();
    });
    ON_CALL(*this, residual).WillByDefault(
      [this](const state_type& x, residual_type& R) {
      mySystem_.residual(x, R);
    });
    ON_CALL(*this, jacobian).WillByDefault(
      [this](const state_type& x, jacobian_type& J) {
      mySystem_.jacobian(x, J);
    });
  }

  MOCK_CONST_METHOD0(createState,    state_type());
  MOCK_CONST_METHOD0(createResidual, residual_type());
  MOCK_CONST_METHOD0(createJacobian, jacobian_type());
  MOCK_CONST_METHOD2(residual, void(const state_type& x, residual_type& R));
  MOCK_CONST_METHOD2(jacobian, void(const state_type& x, jacobian_type& J));

  const concrete_t & viewConcrete(){ return mySystem_; }

 private:
  concrete_t mySystem_;
};

// ==============================================
// problem class with fused residual-jacobian API 
// ==============================================
template<bool use_dense_jacobian = false>
struct ConcreteSystemRJFused
{
  using helper_t = ConcreteSystemRJ<use_dense_jacobian>;
  using scalar_type   = typename helper_t::scalar_type;
  using state_type    = typename helper_t::state_type;
  using residual_type = typename helper_t::residual_type;
  using jacobian_type = typename helper_t::jacobian_type;

 public:
  state_type createState() const { return mySystem_.createState(); }
  residual_type createResidual() const { return mySystem_.createResidual(); }
  jacobian_type createJacobian() const { return mySystem_.createJacobian(); }
  void residualAndJacobian(const state_type& x, residual_type& res, 
    jacobian_type& jac, const bool /*recomputeJacobian*/) const
  {
    mySystem_.residual(x, res);
    mySystem_.jacobian(x, jac);
  }

private:
  helper_t mySystem_;
};

template<bool use_dense_jacobian = false>
struct MockSystemRJFused
{
  using concrete_t = ConcreteSystemRJFused<use_dense_jacobian>;
  using scalar_type   = typename concrete_t::scalar_type;
  using state_type    = typename concrete_t::state_type;
  using residual_type = typename concrete_t::residual_type;
  using jacobian_type = typename concrete_t::jacobian_type;

 public:
  MockSystemRJFused()
  {
    ON_CALL(*this, createState).WillByDefault(
      [this]() {
      return mySystem_.createState();
    });
    ON_CALL(*this, createResidual).WillByDefault(
      [this]() {
      return mySystem_.createResidual();
    });
    ON_CALL(*this, createJacobian).WillByDefault(
      [this]() {
      return mySystem_.createJacobian();
    });
    ON_CALL(*this, residualAndJacobian).WillByDefault(
      [this](const state_type& x, residual_type & R, jacobian_type& J, bool computeJ) {
      mySystem_.residualAndJacobian(x, R, J, computeJ);
    });
  }

  MOCK_CONST_METHOD0(createState, state_type());
  MOCK_CONST_METHOD0(createResidual, residual_type());
  MOCK_CONST_METHOD0(createJacobian, jacobian_type());
  MOCK_CONST_METHOD4(residualAndJacobian, void(const state_type& x, residual_type& R, 
    jacobian_type& J, const bool recomputeJacobian));

  const concrete_t & viewConcrete(){ return mySystem_; }

 private:
  concrete_t mySystem_;
};


// ==============================================
// problem class with hessian-gradient API 
// ==============================================
struct ConcreteSystemHessGrad
{
  using scalar_type = double;
  using state_type  = Eigen::VectorXd;
  using gradient_type = state_type;
  using hessian_type = Eigen::MatrixXd;
  using residual_norm_type = double;

 public:
  state_type createState() const { return state_type(5); }
  gradient_type createGradient() const { return gradient_type(5); }
  hessian_type  createHessian() const  { return hessian_type(5,5); }
  void gradient(const state_type& /*x*/, 
                gradient_type& g, 
                pressio::Norm /*whichNorm*/, 
                residual_norm_type & /*residualNorm*/,
                bool /*recomputeJacobian*/) const
  {
    g.setConstant(3.4);
  }

  void residualNorm(const state_type& /*x*/, 
                pressio::Norm /*whichNorm*/, 
                residual_norm_type & residualNorm) const
  {
    residualNorm = 8.;
  }

  void hessian(const state_type& /*x*/, hessian_type& H) const
  {
    H = goldHessian();
  }

  hessian_type goldHessian() const
  {
    hessian_type H(5,5);
    H.setZero();
    H(0,1) = 1.; H(0,2) = 4.;  H(0,4) = 2.;
    H(1,2) = 3.; H(1,3) = -2.; H(1,4) = 2.;
    H(2,0) = 1.; H(2,1) = 1.;  H(2,4) = 2.;
    H(3,0) = 4.; H(3,1) = -1.; H(3,3) = 1.; H(3,4) = 2.; 
    H(4,0) = 3.;  H(4,3) = 1.; H(4,4) = 5.;
    return H;
  }
};

struct MockSystemHessGrad
{
  using concrete_t = ConcreteSystemHessGrad;
  using scalar_type   = typename concrete_t::scalar_type;
  using state_type    = typename concrete_t::state_type;
  using gradient_type = typename concrete_t::gradient_type;
  using hessian_type = typename concrete_t::hessian_type;
  using residual_norm_type = double;

 public:
  MockSystemHessGrad()
  {
    ON_CALL(*this, createState).WillByDefault(
      [this]() {
      return mySystem_.createState();
    });
    ON_CALL(*this, createGradient).WillByDefault(
      [this]() {
      return mySystem_.createGradient();
    });
    ON_CALL(*this, createHessian).WillByDefault(
      [this]() {
      return mySystem_.createHessian();
    });
    ON_CALL(*this, gradient).WillByDefault(
      [this](const state_type& x, gradient_type& g, pressio::Norm whichNorm, 
        residual_norm_type & residualNorm, bool recomputeJacobian) 
      {
      mySystem_.gradient(x, g, whichNorm, residualNorm, recomputeJacobian);
    });
    ON_CALL(*this, hessian).WillByDefault(
      [this](const state_type& x, hessian_type& H) {
      mySystem_.hessian(x, H);
    });
    ON_CALL(*this, residualNorm).WillByDefault(
      [this](const state_type& x, pressio::Norm whichNorm, residual_norm_type & residualNorm) 
      {
      mySystem_.residualNorm(x, whichNorm, residualNorm);
    });
  }

  MOCK_CONST_METHOD0(createState, state_type());
  MOCK_CONST_METHOD0(createGradient, gradient_type());
  MOCK_CONST_METHOD0(createHessian, hessian_type());
  MOCK_CONST_METHOD5(gradient, void(const state_type&, gradient_type&, pressio::Norm, residual_norm_type&, bool));
  MOCK_CONST_METHOD2(hessian, void(const state_type&, hessian_type&));
  MOCK_CONST_METHOD3(residualNorm, void(const state_type, pressio::Norm, residual_norm_type&));

  const concrete_t & viewConcrete(){ return mySystem_; }

 private:
  concrete_t mySystem_;
};


// ==============================================
// problem class with fused hessian-gradient API 
// ==============================================
struct ConcreteSystemHessGradFused
{
  using scalar_type = double;
  using state_type  = Eigen::VectorXd;
  using gradient_type = state_type;
  using hessian_type = Eigen::MatrixXd;
  using residual_norm_type = double;

 public:
  state_type createState() const { return mySystem_.createState(); }
  gradient_type createGradient() const { return mySystem_.createGradient(); }
  hessian_type  createHessian() const  { return mySystem_.createHessian(); }
  void hessianAndGradient(const state_type& x, 
                hessian_type & H,
                gradient_type& g, 
                pressio::Norm /*whichNorm*/, 
                residual_norm_type & /*residualNorm*/,
                bool recomputeJacobian) const
  {
    g.setConstant(3.4);
    if (recomputeJacobian){
      mySystem_.hessian(x, H);
    }
  }

  void residualNorm(const state_type& /*x*/, 
                pressio::Norm /*whichNorm*/, 
                residual_norm_type & residualNorm) const
  {
    residualNorm = 8.;
  }

  hessian_type goldHessian() const
  {
    return mySystem_.goldHessian();
  }

private:
  ConcreteSystemHessGrad mySystem_;
};

struct MockSystemHessGradFused
{
  using concrete_t = ConcreteSystemHessGradFused;
  using scalar_type   = typename concrete_t::scalar_type;
  using state_type    = typename concrete_t::state_type;
  using gradient_type = typename concrete_t::gradient_type;
  using hessian_type = typename concrete_t::hessian_type;
  using residual_norm_type = double;

 public:
  MockSystemHessGradFused()
  {
    ON_CALL(*this, createState).WillByDefault(
      [this]() {
      return mySystem_.createState();
    });
    ON_CALL(*this, createGradient).WillByDefault(
      [this]() {
      return mySystem_.createGradient();
    });
    ON_CALL(*this, createHessian).WillByDefault(
      [this]() {
      return mySystem_.createHessian();
    });
    ON_CALL(*this, hessianAndGradient).WillByDefault(
      [this](const state_type& x, hessian_type &H, gradient_type& g, pressio::Norm whichNorm, 
        residual_norm_type & residualNorm, bool recomputeJacobian) 
      {
      mySystem_.hessianAndGradient(x, H, g, whichNorm, residualNorm, recomputeJacobian);
    });
    ON_CALL(*this, residualNorm).WillByDefault(
      [this](const state_type& x, pressio::Norm whichNorm, residual_norm_type & residualNorm) 
      {
      mySystem_.residualNorm(x, whichNorm, residualNorm);
    });
  }

  MOCK_CONST_METHOD0(createState, state_type());
  MOCK_CONST_METHOD0(createGradient, gradient_type());
  MOCK_CONST_METHOD0(createHessian, hessian_type());
  MOCK_CONST_METHOD6(hessianAndGradient, 
    void(const state_type&, hessian_type &, gradient_type&, pressio::Norm, residual_norm_type&, bool));
  MOCK_CONST_METHOD3(residualNorm, void(const state_type, pressio::Norm, residual_norm_type&));

  const concrete_t & viewConcrete(){ return mySystem_; }

 private:
  concrete_t mySystem_;
};

// ==============================================
// weighting functor  
// ==============================================
struct ConcreteWeightingOperator
{
  using rank1_operand_type = Eigen::VectorXd;
  using rank2_operand_type = Eigen::MatrixXd;

  void operator()(const rank1_operand_type & /*operand*/, 
                  rank1_operand_type & result) const
  {
    result.setConstant(3.);
  }
  void operator()(const rank2_operand_type & /*operand*/, 
                  rank2_operand_type & result) const
  {
    result.setConstant(2.2);
  }
};

struct MockWeightingOperator
{
  using concrete_t = ConcreteWeightingOperator;
  using rank1_operand_type = typename concrete_t::rank1_operand_type;
  using rank2_operand_type = typename concrete_t::rank2_operand_type;

 public:
  MockWeightingOperator()
  {
    ON_CALL(*this, compute1).WillByDefault(
      [this](const rank1_operand_type & oin, rank1_operand_type & oout) 
      {
        mySystem_(oin, oout);
      });
    ON_CALL(*this, compute2).WillByDefault(
      [this](const rank2_operand_type & oin, rank2_operand_type & oout) 
      {
        mySystem_(oin, oout);
      });
  }

  MOCK_CONST_METHOD2(compute1, 
      void(const rank1_operand_type & oin, rank1_operand_type & oout));
  MOCK_CONST_METHOD2(compute2, 
      void(const rank2_operand_type & oin, rank2_operand_type & oout));

  void operator()(const rank1_operand_type & oin, rank1_operand_type & oout) const
  {
    compute1(oin, oout);
  }

  void operator()(const rank2_operand_type & oin, rank2_operand_type & oout) const
  {
    compute2(oin, oout);
  }

  const concrete_t & viewConcrete(){ return mySystem_; }

 private:
  concrete_t mySystem_;
};


// ==============================================
// linear solver
// ==============================================
template<class A_t, class b_t, class x_t>
struct MockLinearSolver
{
  using matrix_type = A_t;
 public:
  MOCK_METHOD3(solve, void(const A_t &, const b_t &, x_t &));
};

#endif

