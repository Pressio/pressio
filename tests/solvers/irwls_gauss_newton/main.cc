
#include "common.hpp"

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
    if(iterCountR_==1) R.data()->setConstant(1.);
    if(iterCountR_==2) R.data()->setConstant(4.);
    if(iterCountR_==3) R.data()->setConstant(16.);
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

  std::string checkStr = "PASSED";

  {
    using system_t = MySystem;
    system_t sysObj(checkStr);

    MyLinSolverNormalEq linSolverObj(checkStr);

    using state_t = typename system_t::state_type;
    state_t x(numVars);

    using pressio::solvers::nonlinear::experimental::createIRWGaussNewton;
    auto solver = createIRWGaussNewton(sysObj, x, linSolverObj);
    solver.setMaxIterations(3);
    solver.setStoppingCriterion(pressio::solvers::nonlinear::stop::afterMaxIters);
    solver.solve(sysObj, x);
  }

  std::cout << checkStr << std::endl;
  return 0;
}
