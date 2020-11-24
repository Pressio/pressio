
#ifndef PRESSIO_TESTS_SOLVERS_FROZEN_JAC_COMMON_HPP_
#define PRESSIO_TESTS_SOLVERS_FROZEN_JAC_COMMON_HPP_

#include "pressio_solvers.hpp"
#include "types.hpp"

constexpr int numEquations = 5;
constexpr int numVars = 3;

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
    std::cout << *H.data() << std::endl;
    std::cout << *g.data() << std::endl;

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
      if( !H.data()->isApprox(trueH1) ) checkStr_="FAILED";
      if( !g.data()->isApprox(trueG1) ) checkStr_="FAILED";
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
      if( !H.data()->isApprox(trueH2) ) checkStr_="FAILED";
      if( !g.data()->isApprox(trueG2) ) checkStr_="FAILED";
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
      if( !H.data()->isApprox(trueH3) ) checkStr_="FAILED";
      if( !g.data()->isApprox(trueG3) ) checkStr_="FAILED";
    }
  }
};

#endif
