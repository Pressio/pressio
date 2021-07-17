
#include "common.hpp"

struct WeightingOperator
{

#ifdef USE_WRAPPERS
  void operator()(const vec_type & operand, vec_type & result) const{
    result.setConstant(3.);
  }
  void operator()(const mat_type & operand, mat_type & result) const{
    result.setConstant(2.2);
  }

#elif USE_NATIVE
  void operator()(const eig_vec & operand, eig_vec & result) const{
    result.setConstant(3.);
  }
  void operator()(const eig_mat & operand, eig_mat & result) const{
    result.setConstant(2.2);
  }
#endif
};

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

int main()
{
  std::string checkStr = "PASSED";

  {
    using system_t = MySystem;
    system_t sysObj(checkStr);

    MyLinSolverNormalEq linSolverObj(checkStr);

    using state_t = typename system_t::state_type;
    state_t x(numVars);

#if defined USE_GN_NEQ
    auto solver = pressio::nonlinearsolvers::createGaussNewton
      (sysObj,x,linSolverObj, WeightingOperator());
    solver.setMaxIterations(2);
    solver.setStoppingCriterion(pressio::nonlinearsolvers::stop::afterMaxIters);
    solver.solve(sysObj, x);
#endif

#if defined USE_LM_NEQ
    auto solver = pressio::nonlinearsolvers::createLevenbergMarquardt
      (sysObj,x,linSolverObj, WeightingOperator());
    solver.setMaxIterations(2);
    solver.setStoppingCriterion(pressio::nonlinearsolvers::stop::afterMaxIters);
    solver.solve(sysObj, x);
#endif
  }

  std::cout << checkStr << std::endl;
  return 0;
}
