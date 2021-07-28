
#include "pressio_solvers.hpp"
#include "common.hpp"

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
                           bool updateJacobian) const
  {
    ++iterCountR_;
    for (auto i=0; i<R.size(); ++i){
      R_(i) = R(i) + 1.;
    }

    pressio::ops::deep_copy(R, R_);
    if (updateJacobian)
    {
      ++iterCountJ_;
    
      if (iterCountJ_==1){
        jac.setConstant(1.);    
      } 
      else{
        for (auto i=0; i<jac.rows(); ++i)
          for (auto j=0; j<jac.cols(); ++j){
            jac(i,j) += 1.;
          }
      }
    }
  }
};

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
    x.setZero();

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
    auto solver = pressio::nonlinearsolvers::createGaussNewtonQR(sysObj, x, solverObjQR);    
    solver.setMaxIterations(6);
    solver.setStoppingCriterion(pressio::nonlinearsolvers::stop::afterMaxIters);
    solver.setSystemJacobianUpdateFreq(3);
    solver.solve(sysObj, x);

#endif

// test GN normal eq
#if defined USE_GN_NEQ
    auto solver = pressio::nonlinearsolvers::createGaussNewton(sysObj, x, linSolverObj);    
    solver.setMaxIterations(6);
    solver.setStoppingCriterion(pressio::nonlinearsolvers::stop::afterMaxIters);
    solver.setSystemJacobianUpdateFreq(3);
    solver.solve(sysObj, x);
#endif

// test LM 
#if defined USE_LM_NEQ
    auto solver = pressio::nonlinearsolvers::createLevenbergMarquardt(sysObj, x, linSolverObj);    
    solver.setMaxIterations(6);
    solver.setStoppingCriterion(pressio::nonlinearsolvers::stop::afterMaxIters);
    solver.setSystemJacobianUpdateFreq(3);
    solver.solve(sysObj, x);
#endif
  }

  std::cout << checkStr << std::endl;
  return 0;
}
