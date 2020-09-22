
#include "pressio_solvers.hpp"
#include "common.hpp"

struct MySystem
{
  using scalar_type = double;
  using state_type = vec_type;
  using hessian_type  = mat_type;
  using gradient_type = vec_type;
  using residual_type = state_type;
  using jacobian_type = mat_type;

  std::string & checkStr_;
  mutable std::size_t iterCountR_ = {};
  mutable std::size_t iterCountJ_ = {};
  mutable residual_type R_;
  mutable jacobian_type J_;

  MySystem(std::string & checkStr) 
    : checkStr_(checkStr), 
      R_(numEquations),
      J_(numEquations, numVars){}

  ~MySystem()
  {
    // when this is destructed, 
    // we know what iterCount_ should be 
    if(iterCountR_ != 6) checkStr_="FAILED";
    if(iterCountJ_ != 6) checkStr_="FAILED";
  }

  hessian_type createHessian() const {
    return hessian_type(numVars, numVars);
  }

  gradient_type createGradient() const {
    return gradient_type(numVars);
  }

  void hessian(const state_type &, 
               hessian_type & H) const
  {
    computeJacobian();

    *H.data() =  J_.data()->transpose() * (*J_.data());
  }

  void gradient(const state_type & , 
                gradient_type & g, 
                pressio::Norm normKind, 
                scalar_type & residualNorm,
                bool updateJacobian) const
  {
    ++iterCountR_;

    for (auto i=0; i<R_.extent(0); ++i)
      R_[i] += 1.;

    if (normKind == pressio::Norm::L2)
      residualNorm = R_.data()->norm();

    if (updateJacobian){
      computeJacobian();
    }

    *g.data() = J_.data()->transpose() *  (*R_.data());
  }

  void residualNorm(const state_type &, 
                    pressio::Norm normKind, 
                    scalar_type & resNorm) const
  {
    if (normKind == pressio::Norm::L2)
      resNorm = R_.data()->norm();
  }

private:
  void computeJacobian() const
  {  
    ++iterCountJ_;
    if (iterCountJ_==1 or iterCountJ_==2){
      J_.data()->setConstant(1.);
    } 
    if (iterCountJ_==3  or iterCountJ_==5){
      for (auto i=0; i<J_.extent(0); ++i)
        for (auto j=0; j<J_.extent(1); ++j)
          J_(i,j) += 1.;
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

// test GN normal eq
#if defined USE_GN_NEQ
    using solver = pressio::solvers::nonlinear::composeGaussNewton_t<
      system_t, 
      // pressio::solvers::nonlinear::DefaultUpdate,
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
      system_t, 
      // pressio::solvers::nonlinear::DefaultUpdate,
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
