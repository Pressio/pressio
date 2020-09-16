
#include "pressio_solvers.hpp"
#include "common.hpp"

// using eig_mat = Eigen::Matrix<double, -1, -1>;
// using eig_vec = Eigen::VectorXd;
// using vec_type = pressio::containers::Vector<eig_vec>;
// using mat_type = pressio::containers::Matrix<eig_mat>;
// constexpr int numEquations = 5;
// constexpr int numVars = 3;

struct MySystem
{
  using scalar_type = double;
  using state_type = vec_type;
  using residual_type = state_type;
  using jacobian_type = mat_type;

private:
  std::string & checkStr_;
  mutable std::size_t iterCountR_ = {};
  mutable std::size_t iterCountJ_ = {};
  mutable residual_type R_;

public:
  MySystem(std::string & checkStr) 
    : checkStr_(checkStr), 
      R_(numEquations){}

  ~MySystem()
  {
    // when this is destructed, 
    // we know what iterCount_ should be 
    if(iterCountR_ != 6) checkStr_="FAILED";
    if(iterCountJ_ != 3) checkStr_="FAILED";
  }

  residual_type createResidual() const {
    return residual_type(numEquations);
  }

  jacobian_type createJacobian() const {
    return jacobian_type(numEquations, numVars);
  }

  void residualAndJacobian(const state_type & y,
                           residual_type & R, 
                           jacobian_type & jac,
                           pressio::Norm normKind, 
                           scalar_type & resNorm, 
                           bool updateJacobian) const
  {
    ++iterCountR_;
    for (auto i=0; i<R.extent(0); ++i)
      R_[i] = R[i] + 1.;

    pressio::ops::deep_copy(R, R_);
    if (normKind == pressio::Norm::L2)
      resNorm = R_.data()->norm();

    if (updateJacobian)
    {
      ++iterCountJ_;
    
      if (iterCountJ_==1){
        jac.data()->setConstant(1.);    
      } 
      else{
        for (auto i=0; i<jac.extent(0); ++i)
          for (auto j=0; j<jac.extent(1); ++j)
            jac(i,j) += 1.;
      }
    }
  }

  void residualNorm(const state_type &, 
                    pressio::Norm normKind, 
                    scalar_type & resNorm) const
  {
    if (normKind == pressio::Norm::L2)
      resNorm = R_.data()->norm();
  }
};

// struct MyLinSolverNormalEq
//   : pressio::solvers::LinearBase<mat_type, MyLinSolverNormalEq>
// {
//   using matrix_type = mat_type;

//   std::size_t iterCount_ = {};
//   std::string & checkStr_;

//   MyLinSolverNormalEq(std::string & checkStr) 
//   : checkStr_(checkStr){}

//   template<typename H_t, typename state_t>
//   void solve(const H_t & H, 
//              const state_t & g, 
//              state_t & correction)
//   {
//     ++iterCount_;
//     std::cout << iterCount_ << std::endl;
//     std::cout << *H.data() << std::endl;
//     std::cout << *g.data() << std::endl;

//     eig_mat trueH1(numVars,numVars);
//     trueH1 << 5., 5., 5., 5., 5., 5., 5., 5., 5.;
// #if defined USE_LM_NEQ
//     // LM adds mu*diag(H)
//     trueH1(0,0) += 5.; trueH1(1,1) += 5.; trueH1(2,2) += 5.; 
// #endif

//     eig_vec trueG1(numVars);
//     trueG1 << -5., -5., -5.;

//     eig_mat trueH2 = trueH1;
//     eig_vec trueG2(numVars);
//     trueG2 << -10., -10., -10.;

//     eig_mat trueH3(numVars,numVars);
//     trueH3 << 20., 20., 20., 20., 20., 20., 20., 20., 20.;
// #if defined USE_LM_NEQ
//     // LM adds mu*diag(H)
//     trueH3(0,0) += 20.; trueH3(1,1) += 20.; trueH3(2,2) += 20.; 
// #endif
//     eig_vec trueG3(numVars);
//     trueG3 << -30., -30., -30.;

//     eig_mat trueH4 = trueH3;
//     eig_vec trueG4(numVars);
//     trueG4 << -40., -40., -40.;

//     eig_mat trueH5 = trueH3;
//     eig_vec trueG5(numVars);
//     trueG5 << -50., -50., -50.;

//     eig_mat trueH6(numVars,numVars);
//     trueH6 << 45., 45., 45., 45., 45., 45., 45., 45., 45.;
// #if defined USE_LM_NEQ
//     // LM adds mu*diag(H)
//     trueH6(0,0) += 45.; trueH6(1,1) += 45.; trueH6(2,2) += 45.; 
// #endif
//     eig_vec trueG6(numVars);
//     trueG6 << -90., -90., -90.;

//     if (iterCount_==1){ 
//       if( ! H.data()->isApprox(trueH1) ) checkStr_="FAILED";
//       if( ! g.data()->isApprox(trueG1) ) checkStr_="FAILED";
//     }

//     if (iterCount_==2){ 
//       if( ! H.data()->isApprox(trueH2) ) checkStr_="FAILED";
//       if( ! g.data()->isApprox(trueG2) ) checkStr_="FAILED";
//     }

//     if (iterCount_==3){ 
//       if( ! H.data()->isApprox(trueH3) ) checkStr_="FAILED";
//       if( ! g.data()->isApprox(trueG3) ) checkStr_="FAILED";
//     }

//     if (iterCount_==4){ 
//       if( ! H.data()->isApprox(trueH4) ) checkStr_="FAILED";
//       if( ! g.data()->isApprox(trueG4) ) checkStr_="FAILED";
//     }

//     if (iterCount_==5){ 
//       if( ! H.data()->isApprox(trueH5) ) checkStr_="FAILED";
//       if( ! g.data()->isApprox(trueG5) ) checkStr_="FAILED";
//     }

//     if (iterCount_==6){ 
//       if( ! H.data()->isApprox(trueH6) ) checkStr_="FAILED";
//       if( ! g.data()->isApprox(trueG6) ) checkStr_="FAILED";
//     }
//   }  
// };

// struct MySolverQR
// {
//   std::size_t iterCount_ = {};
//   std::string & checkStr_;

//   MySolverQR(std::string & checkStr) 
//   : checkStr_(checkStr){}

//   template<typename r_t, typename result_t>
//   void applyQTranspose(const r_t & r, result_t & QTr)
//   {}

//   template<typename result_t, typename g_t>
//   void applyRTranspose(result_t & QTr, g_t & g)
//   {}

//   template<typename qt_t, typename c_t>
//   void solve(const qt_t & QTr, c_t & correction)
//   {}

//   template<typename J_t>
//   void computeThin(const J_t & Jac)
//   {
//     ++iterCount_;
//     std::cout << Jac.extent(0) << " " << Jac.extent(1) << std::endl;
//     // std::cout << *H.data() << std::endl;
//     // std::cout << *g.data() << std::endl;

//     eig_mat trueJ1(numEquations,numVars);
//     trueJ1.setConstant(1.);

//     eig_mat trueJ2 = trueJ1;

//     eig_mat trueJ3(numEquations,numVars);
//     trueJ3.setConstant(2.);

//     eig_mat trueJ4 = trueJ3;

//     eig_mat trueJ5 = trueJ3;

//     eig_mat trueJ6(numEquations,numVars);
//     trueJ6.setConstant(3.);

//     if (iterCount_==1){ 
//       if( ! Jac.data()->isApprox(trueJ1) ) checkStr_="FAILED";
//     }

//     if (iterCount_==2){ 
//       if( ! Jac.data()->isApprox(trueJ2) ) checkStr_="FAILED";
//     }

//     if (iterCount_==3){ 
//       if( ! Jac.data()->isApprox(trueJ3) ) checkStr_="FAILED";
//     }

//     if (iterCount_==4){ 
//       if( ! Jac.data()->isApprox(trueJ4) ) checkStr_="FAILED";
//     }

//     if (iterCount_==5){ 
//       if( ! Jac.data()->isApprox(trueJ5) ) checkStr_="FAILED";
//     }

//     if (iterCount_==6){ 
//       if( ! Jac.data()->isApprox(trueJ6) ) checkStr_="FAILED";
//     }
//   }  
// };

// namespace pressio{ namespace qr{  namespace details{ 

// template<>
// struct traits<MySolverQR>{
//   using matrix_t = mat_type;  
//   using concrete_t = double;  
// };
// }}}


int main()
{
  std::string checkStr = "PASSED";

  {
    using system_t = MySystem;
    system_t sysObj(checkStr);

#if defined USE_GN_NEQ or USE_LM_NEQ
    MyLinSolverNormalEq linSolverObj(checkStr);
#endif
#if defined USE_GN_QR
    MySolverQR solverObjQR(checkStr);
#endif

    using state_t = typename system_t::state_type;
    state_t x(numVars);

    /*
    what we want to test here? 
    We want to test that the frozen Jacobian behaves correctly. 
    To do this, we don't solve a real problem but we just want 
    to test tha the update of jacobin is done correctly for  
    any of the methods possible.
    We test this by keeping track of the calls to the 
    solver to make sure the operators are what we expect 
    for given iterations.

    - we fix num of iterations to 6 
    - update jacobian every 3 steps, 
      so it should be updated at step 3 and 6

    - we start from J[:,:]=1, 
      at first update we have J[:,:]=2

    - residual starts from R[:]=1, incremented 
      element-wise by 1 at every call

    - with these assumptions, we should know what the 
      operators look like at every call
     */

// test GN with QR
#if defined USE_GN_QR
    using solver = pressio::solvers::nonlinear::composeGaussNewtonQR_t<
      system_t, pressio::solvers::nonlinear::DefaultUpdate,
      MySolverQR>;
    solver solver1(sysObj, x, solverObjQR);
    solver1.setMaxIterations(6);
    solver1.setStoppingCriterion(pressio::solvers::nonlinear::stop::afterMaxIters);
    solver1.setSystemJacobianUpdateFreq(3);
    solver1.solve(sysObj, x);

#endif

// test GN normal eq
#if defined USE_GN_NEQ
    using solver = pressio::solvers::nonlinear::composeGaussNewton_t<
      system_t, pressio::solvers::nonlinear::DefaultUpdate,
      MyLinSolverNormalEq>;
    solver solver1(sysObj, x, linSolverObj);
    solver1.setMaxIterations(6);
    solver1.setStoppingCriterion(pressio::solvers::nonlinear::stop::afterMaxIters);
    solver1.setSystemJacobianUpdateFreq(3);
    solver1.solve(sysObj, x);
#endif

// test LM 
#if defined USE_LM_NEQ
    using solver = pressio::solvers::nonlinear::composeLevenbergMarquardt_t<
      system_t, pressio::solvers::nonlinear::DefaultUpdate,
      MyLinSolverNormalEq>;
    solver solver1(sysObj, x, linSolverObj);
    solver1.setMaxIterations(6);
    solver1.setStoppingCriterion(pressio::solvers::nonlinear::stop::afterMaxIters);
    solver1.setSystemJacobianUpdateFreq(3);
    solver1.solve(sysObj, x);
#endif
  }

  std::cout << checkStr << std::endl;
  return 0;
}
