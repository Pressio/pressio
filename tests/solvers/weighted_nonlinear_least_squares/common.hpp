
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
    : checkStr_(checkStr) {}

  template <typename H_t, typename state_t>
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
    eig_mat trueH1(numVars, numVars);
    trueH1 << 11., 11., 11., 11., 11., 11., 11., 11., 11.;
#if defined USE_LM_NEQ
    // LM adds mu*diag(H)
    trueH1(0, 0) += 11.;
    trueH1(1, 1) += 11.;
    trueH1(2, 2) += 11.;
#endif

    // G1 = J1^T Mr
    // where J1[:] = 1
    // where Mr[:] = 3
    eig_vec trueG1(numVars);
    trueG1 << -15., -15., -15.;

    if(iterCount_ == 1) {
      if(!H.data()->isApprox(trueH1))
	checkStr_ = "FAILED";
      if(!g.data()->isApprox(trueG1))
	checkStr_ = "FAILED";
    }


    ///////////////////////////
    /////// second call  ///////
    ///////////////////////////
    // H2 = J2 * MJ2
    // where J2[:] = 2
    // where MJ2[:] = 2.2
    eig_mat trueH2(numVars, numVars);
    trueH2 << 22., 22., 22., 22., 22., 22., 22., 22., 22.;
#if defined USE_LM_NEQ
    // LM adds mu*diag(H)
    trueH2(0, 0) += 22.;
    trueH2(1, 1) += 22.;
    trueH2(2, 2) += 22.;
#endif

    // G2 = J2^T Mr2
    // where J2[:] = 2
    // where Mr2[:] = 3
    eig_vec trueG2(numVars);
    trueG2 << -30., -30., -30.;

    if(iterCount_ == 2) {
      if(!H.data()->isApprox(trueH2))
	checkStr_ = "FAILED";
      if(!g.data()->isApprox(trueG2))
	checkStr_ = "FAILED";
    }
  }
};

#endif
