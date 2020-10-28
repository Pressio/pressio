
#include "pressio_solvers.hpp"
#include "common.hpp"

struct WeightingOperator
{
  using residual_type = vec_type;
  using jacobian_type = mat_type;

  void operator()(const residual_type & operand,
                  residual_type & result) const
  {
    // fake the operation
    result.data()->setConstant(3.);
  }

  void operator()(const jacobian_type & operand,
                  jacobian_type & result) const
  {
    // fake the operation
    result.data()->setConstant(2.2);
  }
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
    a.data()->setConstant(0);
    return a;
  }

  jacobian_type createJacobian() const {
    auto a = jacobian_type(numEquations, numVars);
    a.data()->setConstant(0);
    return a;
  }

  void residual(const state_type& x,
    residual_type & R) const
  {
    ++iterCountR_;

    for (auto i=0; i<R.extent(0); ++i)
      R[i] += 1.;

    // if (normKind == pressio::Norm::L2)
    //   normResidual = R.data()->norm();
  }

  void jacobian(const state_type& x, jacobian_type & jac) const
  {
    ++iterCountJ_;

    for (auto i=0; i<jac.extent(0); ++i)
      for (auto j=0; j<jac.extent(1); ++j)
        jac(i,j) += 1.;
  }
};

int main()
{
  std::string checkStr = "PASSED";

  {
    using system_t = MySystem;
    system_t sysObj(checkStr);

    MyLinSolverNormalEq linSolverObj(checkStr);
    WeightingOperator M;

    using state_t = typename system_t::state_type;
    state_t x(numVars);

#if defined USE_GN_NEQ
    auto solver = pressio::solvers::nonlinear::createGaussNewton(sysObj,x,linSolverObj,M);
    solver.setMaxIterations(2);
    solver.setStoppingCriterion(pressio::solvers::nonlinear::stop::afterMaxIters);
    solver.solve(sysObj, x);
#endif

#if defined USE_LM_NEQ
    auto solver = pressio::solvers::nonlinear::createLevenbergMarquardt(sysObj,x,linSolverObj,M);
    solver.setMaxIterations(2);
    solver.setStoppingCriterion(pressio::solvers::nonlinear::stop::afterMaxIters);
    solver.solve(sysObj, x);
#endif
  }

  std::cout << checkStr << std::endl;
  return 0;
}
