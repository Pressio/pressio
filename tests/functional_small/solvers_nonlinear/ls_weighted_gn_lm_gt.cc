
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "pressio/solvers.hpp"

using eig_mat = Eigen::Matrix<double, -1, -1>;
using eig_vec = Eigen::VectorXd;

constexpr int numEquations = 5;
constexpr int numVars = 3;

struct MyLinSolverNormalEq
{
  using matrix_type = eig_mat;

  std::size_t iterCount_ = {};
  std::string & checkStr_;
  std::string methodName_;

  MyLinSolverNormalEq(std::string & checkStr, 
                      std::string methodName)
  : checkStr_(checkStr), methodName_(methodName){}

  template<typename H_t, typename state_t>
  void solve(const H_t & H,
             const state_t & g,
             state_t & correction)
  {
    ++iterCount_;
    std::cout << iterCount_ << std::endl;
    std::cout << *H.data() << std::endl;
    std::cout << *g.data() << std::endl;

    ///////////////////////////
    /////// first call  ///////
    ///////////////////////////
    // H1 = J1 * MJ
    // where J1[:] = 1
    // where MJ[:] = 2.2
    eig_mat trueH1(numVars,numVars);
    trueH1 << 11., 11., 11., 11., 11., 11., 11., 11., 11.;
    if (methodName_=="lm"){
      // LM adds mu*diag(H)
      trueH1(0,0) += 11.; trueH1(1,1) += 11.; trueH1(2,2) += 11.;
    }

    // G1 = J1^T Mr
    // where J1[:] = 1
    // where Mr[:] = 3
    eig_vec trueG1(numVars);
    trueG1 << -15., -15., -15.;

    if (iterCount_==1){
      if( ! H.isApprox(trueH1) ) checkStr_="FAILED";
      if( ! g.isApprox(trueG1) ) checkStr_="FAILED";
    }


    ///////////////////////////
    /////// second call  ///////
    ///////////////////////////
    // H2 = J2 * MJ2
    // where J2[:] = 2
    // where MJ2[:] = 2.2
    eig_mat trueH2(numVars,numVars);
    trueH2 << 22., 22., 22., 22., 22., 22., 22., 22., 22.;
    if (methodName_=="lm"){
      // LM adds mu*diag(H)
      trueH2(0,0) += 22.; trueH2(1,1) += 22.; trueH2(2,2) += 22.;
    }

    // G2 = J2^T Mr2
    // where J2[:] = 2
    // where Mr2[:] = 3
    eig_vec trueG2(numVars);
    trueG2 << -30., -30., -30.;

    if (iterCount_==2){
      if( ! H.isApprox(trueH2) ) checkStr_="FAILED";
      if( ! g.isApprox(trueG2) ) checkStr_="FAILED";
    }
  }
};


struct WeightingOperator
{
  void operator()(const eig_vec & operand, eig_vec & result) const{
    result.setConstant(3.);
  }
  void operator()(const eig_mat & operand, eig_mat & result) const{
    result.setConstant(2.2);
  }
};

struct MySystem
{
  using state_type = eig_vec;
  using residual_type = state_type;
  using jacobian_type = eig_mat;

  std::string & checkStr_;
  mutable std::size_t iterCountR_ = {};
  mutable std::size_t iterCountJ_ = {};

  MySystem(std::string & checkStr)
    : checkStr_(checkStr){}

  ~MySystem(){}

  state_type createState() const {
    auto a = state_type(numVars);
    a.setConstant(0);
    return a;
  }

  residual_type createResidual() const {
    auto a = residual_type(numEquations);
    a.setConstant(0);
    return a;
  }

  jacobian_type createJacobian() const {
    auto a = jacobian_type(numEquations, numVars);
    a.setConstant(0);
    return a;
  }

  void residual(const state_type& x,
    residual_type & R) const
  {
    ++iterCountR_;

    for (auto i=0; i<R.size(); ++i){
      R(i) += 1.;
    }
  }

  void jacobian(const state_type& x, jacobian_type & jac) const
  {
    ++iterCountJ_;

    for (auto i=0; i<jac.rows(); ++i)
      for (auto j=0; j<jac.cols(); ++j){
        jac(i,j) += 1.;
      }
  }
};

TEST(solvers_nonlinear, weighted_least_squares_gauss_newton)
{
  std::string checkStr = "PASSED";
  using system_t = MySystem;
  system_t sysObj(checkStr);

  MyLinSolverNormalEq linSolverObj(checkStr, "gn");

  using state_t = typename system_t::state_type;
  state_t x(numVars);

  auto solver = pressio::nonlinearsolvers::create_gauss_newton
    (sysObj, linSolverObj, WeightingOperator());
  solver.setMaxIterations(2);
  solver.setStoppingCriterion(pressio::nonlinearsolvers::Stop::AfterMaxIters);
  solver.solve(sysObj, x);
  ASSERT_TRUE(checkStr == "PASSED");
}

TEST(solvers_nonlinear, weighted_least_squares_lm)
{
  std::string checkStr = "PASSED";
  using system_t = MySystem;
  system_t sysObj(checkStr);

  MyLinSolverNormalEq linSolverObj(checkStr, "lm");

  using state_t = typename system_t::state_type;
  state_t x(numVars);

  auto solver = pressio::nonlinearsolvers::create_levenberg_marquardt
    (sysObj, linSolverObj, WeightingOperator());
  solver.setMaxIterations(2);
  solver.setStoppingCriterion(pressio::nonlinearsolvers::Stop::AfterMaxIters);
  solver.solve(sysObj, x);
  ASSERT_TRUE(checkStr == "PASSED");
}
