
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "pressio/solvers.hpp"

  /*
    this test is for IRWLS p=1 with Gauss-Newton

    numEquations = 5;
    numVars = 3;

    r and J are set by the system just for testing
    purposes, there is no numerical meaning in all this.
    we just make sure the steps are correct.

    *************************************
    ***first iter of nonlinear solver***
    *************************************
    r = [1 1 1 1 1]
    J = [1 1 1;
         1 1 1;
   1 1 1;
   1 1 1;
   1 1 1]

    H = J^T W J   where W=[1 ... 1]
      =[5 5 5;
  5 5 5;
  5 5 5]

    b = -J^T W r
      = -5 -5 -5

    *************************************
    ***second iter of nonlinear solver***
    *************************************
    r = [4 4 4 4 4]
    J = [2 2 2;
         2 2 2;
   2 2 2;
   2 2 2;
   2 2 2]

    now W changes because we are doing second iteration
    W = r^(p-2)
      = 0.25 0.25 0.25 0.25 0.25

    H = J^T W J
      =[5 5 5;
  5 5 5;
  5 5 5]

    b = -J^T W r
      = -10 -10 -10

    *************************************
    ***third iter of nonlinear solver***
    *************************************
    r = [16 16 16 16 16]
    J = [3 3 3;
         3 3 3;
   3 3 3;
   3 3 3;
   3 3 3]

    now W changes again we are doing third iteration
    W = r^(p-2)
      = 0.0625 0.0625 0.0625 0.0625 0.0625

    H = J^T W J
      =[2.8125 2.8125 2.8125;
  2.8125 2.8125 2.8125;
  2.8125 2.8125 2.8125]

    b = -J^T W r
      = -15 -15 -15
   */


constexpr int numEquations = 5;
constexpr int numVars = 3;

using eig_mat = Eigen::Matrix<double, -1, -1>;
using eig_vec = Eigen::VectorXd;
using vec_type = eig_vec;
using mat_type = eig_mat;

struct MySystem
{
  using scalar_type = double;
  using state_type = vec_type;
  using residual_type = state_type;
  using jacobian_type = mat_type;

  std::string & checkStr_;
  mutable std::size_t iterCountR_ = {};
  mutable std::size_t iterCountJ_ = {};

  MySystem(std::string & checkStr)
    : checkStr_(checkStr){}

  ~MySystem(){}

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
    if(iterCountR_==1) R.setConstant(1.);
    if(iterCountR_==2) R.setConstant(4.);
    if(iterCountR_==3) R.setConstant(16.);
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

struct MyLinSolverNormalEq
{
  using matrix_type = mat_type;

  std::size_t iterCount_ = {};
  std::string & checkStr_;

  MyLinSolverNormalEq(std::string & checkStr)
  : checkStr_(checkStr){}

  template<typename H_t, typename state_t>
  void solve(const H_t & H,
             const state_t & g,
             state_t & correction)
  {
    ++iterCount_;
    std::cout << iterCount_ << std::endl;
    std::cout << H << std::endl;
    std::cout << g << std::endl;

    ///////////////////////////
    /////// first call  ///////
    ///////////////////////////
    eig_mat trueH1(numVars,numVars); trueH1.setConstant(5.);
// #if defined USE_LM_NEQ
//     // LM adds mu*diag(H)
//     trueH1(0,0) += 11.; trueH1(1,1) += 11.; trueH1(2,2) += 11.;
// #endif
    eig_vec trueG1(numVars); trueG1.setConstant(-5.);
    if (iterCount_==1){
      if( !H.isApprox(trueH1) ) checkStr_="FAILED";
      if( !g.isApprox(trueG1) ) checkStr_="FAILED";
    }

    ///////////////////////////
    /////// second call  ///////
    ///////////////////////////
    eig_mat trueH2(numVars,numVars); trueH2.setConstant(5.);
// #if defined USE_LM_NEQ
//     // LM adds mu*diag(H)
//     trueH1(0,0) += 11.; trueH1(1,1) += 11.; trueH1(2,2) += 11.;
// #endif
    eig_vec trueG2(numVars); trueG2.setConstant(-10.);
    if (iterCount_==2){
      if( !H.isApprox(trueH2) ) checkStr_="FAILED";
      if( !g.isApprox(trueG2) ) checkStr_="FAILED";
    }

    ///////////////////////////
    /////// third call  ///////
    ///////////////////////////
    eig_mat trueH3(numVars,numVars); trueH3.setConstant(2.8125);
// #if defined USE_LM_NEQ
//     // LM adds mu*diag(H)
//     trueH1(0,0) += 11.; trueH1(1,1) += 11.; trueH1(2,2) += 11.;
// #endif
    eig_vec trueG3(numVars); trueG3.setConstant(-15.);
    if (iterCount_==3){
      if( !H.isApprox(trueH3) ) checkStr_="FAILED";
      if( !g.isApprox(trueG3) ) checkStr_="FAILED";
    }
  }
};

TEST(solvers_nonlinear, irwls_gauss_newton)
{
  std::string checkStr = "PASSED";

  using system_t = MySystem;
  system_t sysObj(checkStr);

  MyLinSolverNormalEq linSolverObj(checkStr);

  using state_t = typename system_t::state_type;
  state_t x(numVars);

  using pressio::nonlinearsolvers::experimental::create_irlw_gauss_newton;
  auto solver = create_irlw_gauss_newton(sysObj, x, linSolverObj);
  solver.setMaxIterations(3);
  solver.setStoppingCriterion(pressio::nonlinearsolvers::Stop::afterMaxIters);
  solver.solve(sysObj, x);

}
