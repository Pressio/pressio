
#ifndef PRESSIO_TESTS_SOLVERS_FROZEN_JAC_COMMON_HPP_
#define PRESSIO_TESTS_SOLVERS_FROZEN_JAC_COMMON_HPP_

#include "pressio_solvers.hpp"

using eig_mat = Eigen::Matrix<double, -1, -1>;
using eig_vec = Eigen::VectorXd;
using vec_type = pressio::containers::Vector<eig_vec>;
using mat_type = pressio::containers::DenseMatrix<eig_mat>;

constexpr int numEquations = 5;
constexpr int numVars = 3;


struct MyLinSolverNormalEq
  // : pressio::solvers::LinearBase<mat_type, MyLinSolverNormalEq>
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

    eig_mat trueH1(numVars,numVars);
    trueH1 << 5., 5., 5., 5., 5., 5., 5., 5., 5.;
#if defined USE_LM_NEQ
    // LM adds mu*diag(H)
    trueH1(0,0) += 5.; trueH1(1,1) += 5.; trueH1(2,2) += 5.; 
#endif

    eig_vec trueG1(numVars);
    trueG1 << -5., -5., -5.;

    eig_mat trueH2 = trueH1;
    eig_vec trueG2(numVars);
    trueG2 << -10., -10., -10.;

    eig_mat trueH3(numVars,numVars);
    trueH3 << 20., 20., 20., 20., 20., 20., 20., 20., 20.;
#if defined USE_LM_NEQ
    // LM adds mu*diag(H)
    trueH3(0,0) += 20.; trueH3(1,1) += 20.; trueH3(2,2) += 20.; 
#endif
    eig_vec trueG3(numVars);
    trueG3 << -30., -30., -30.;

    eig_mat trueH4 = trueH3;
    eig_vec trueG4(numVars);
    trueG4 << -40., -40., -40.;

    eig_mat trueH5 = trueH3;
    eig_vec trueG5(numVars);
    trueG5 << -50., -50., -50.;

    eig_mat trueH6(numVars,numVars);
    trueH6 << 45., 45., 45., 45., 45., 45., 45., 45., 45.;
#if defined USE_LM_NEQ
    // LM adds mu*diag(H)
    trueH6(0,0) += 45.; trueH6(1,1) += 45.; trueH6(2,2) += 45.; 
#endif
    eig_vec trueG6(numVars);
    trueG6 << -90., -90., -90.;

    if (iterCount_==1){ 
      if( ! H.data()->isApprox(trueH1) ) checkStr_="FAILED";
      if( ! g.data()->isApprox(trueG1) ) checkStr_="FAILED";
    }

    if (iterCount_==2){ 
      if( ! H.data()->isApprox(trueH2) ) checkStr_="FAILED";
      if( ! g.data()->isApprox(trueG2) ) checkStr_="FAILED";
    }

    if (iterCount_==3){ 
      if( ! H.data()->isApprox(trueH3) ) checkStr_="FAILED";
      if( ! g.data()->isApprox(trueG3) ) checkStr_="FAILED";
    }

    if (iterCount_==4){ 
      if( ! H.data()->isApprox(trueH4) ) checkStr_="FAILED";
      if( ! g.data()->isApprox(trueG4) ) checkStr_="FAILED";
    }

    if (iterCount_==5){ 
      if( ! H.data()->isApprox(trueH5) ) checkStr_="FAILED";
      if( ! g.data()->isApprox(trueG5) ) checkStr_="FAILED";
    }

    if (iterCount_==6){ 
      if( ! H.data()->isApprox(trueH6) ) checkStr_="FAILED";
      if( ! g.data()->isApprox(trueG6) ) checkStr_="FAILED";
    }
  }  
};

 
struct MySolverQR
{
  std::size_t iterCount_ = {};
  std::string & checkStr_;

  MySolverQR(std::string & checkStr) 
  : checkStr_(checkStr){}

  template<typename r_t, typename result_t>
  void applyQTranspose(const r_t & r, result_t & QTr) const
  {}

  template<typename result_t, typename g_t>
  void applyRTranspose(result_t & QTr, g_t & g) const
  {}

  template<typename qt_t, typename c_t>
  void solve(const qt_t & QTr, c_t & correction) const
  {}

  template<typename J_t>
  void computeThin(const J_t & Jac)
  {
    ++iterCount_;
    std::cout << Jac.extent(0) << " " << Jac.extent(1) << std::endl;
    // std::cout << *H.data() << std::endl;
    // std::cout << *g.data() << std::endl;

    eig_mat trueJ1(numEquations,numVars);
    trueJ1.setConstant(1.);

    eig_mat trueJ2 = trueJ1;

    eig_mat trueJ3(numEquations,numVars);
    trueJ3.setConstant(2.);

    eig_mat trueJ4 = trueJ3;

    eig_mat trueJ5 = trueJ3;

    eig_mat trueJ6(numEquations,numVars);
    trueJ6.setConstant(3.);

    if (iterCount_==1){ 
      if( ! Jac.data()->isApprox(trueJ1) ) checkStr_="FAILED";
    }

    if (iterCount_==2){ 
      if( ! Jac.data()->isApprox(trueJ2) ) checkStr_="FAILED";
    }

    if (iterCount_==3){ 
      if( ! Jac.data()->isApprox(trueJ3) ) checkStr_="FAILED";
    }

    if (iterCount_==4){ 
      if( ! Jac.data()->isApprox(trueJ4) ) checkStr_="FAILED";
    }

    if (iterCount_==5){ 
      if( ! Jac.data()->isApprox(trueJ5) ) checkStr_="FAILED";
    }

    if (iterCount_==6){ 
      if( ! Jac.data()->isApprox(trueJ6) ) checkStr_="FAILED";
    }
  }  
};

namespace pressio{ namespace qr{  namespace details{ 

template<>
struct traits<MySolverQR>{
  using matrix_t = mat_type;  
  using concrete_t = double;  
};
}}}

#endif

