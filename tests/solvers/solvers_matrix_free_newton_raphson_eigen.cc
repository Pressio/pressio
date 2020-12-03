
#include "pressio_solvers.hpp"

using eig_v_t = Eigen::VectorXd;

struct ValidSystem
{
private:
  using matrix_n_t = Eigen::SparseMatrix<double>;
  using jacobian_type = pressio::containers::SparseMatrix<matrix_n_t>;

public:
  using scalar_type   = double;
  using state_type    = pressio::containers::Vector<eig_v_t>;
  using residual_type = state_type;

  residual_type createResidual() const {
    return residual_type(2);
  }

  void residual(const state_type & x,
                residual_type & res) const
  {
    res(0) =  x(0)*2.;
    res(1) =  x(1)*2.;
  }

  void applyJacobian(const state_type & x,
		     const residual_type & in,
		     residual_type & out,
		     bool updateJacobianAction) const
  {
    // in theory here we should compute action without
    // computing the jacobian, but here we jsut test things
    jacobian_type jac(2,2);
    jac.data()->coeffRef(0, 0) = x(0);
    jac.data()->coeffRef(0, 1) = x(1);
    jac.data()->coeffRef(1, 0) = x(0)*2.;
    jac.data()->coeffRef(1, 1) = x(1)*2.;
    *out.data() = (*jac.data()) * (*in.data());
  }
};

struct mylinsolver
{
  using state_type    = pressio::containers::Vector<eig_v_t>;
  using residual_type = state_type;

  std::size_t iterCount_ = {};
  std::string & checkStr_;

  mylinsolver(std::string & checkStr)
  : checkStr_(checkStr){}

  // solve Ax = r
  // where A is an object that evaluates the action of the jacobian
  // keep in mind that for this test, this linear solve is called
  // at each newton iteration so x is effectively the correction
  template<class T>
  void solve(const T & jAction,
             const residual_type & r,
             state_type & x)
  {
    ++iterCount_;
    std::cout << iterCount_ << std::endl;
    std::cout << *x.data() << std::endl;
    // std::cout << *r.data() << std::endl;

    eig_v_t trueX(2);
    eig_v_t trueR(2);
    eig_v_t trueJb(2);
    eig_v_t trueJr(2);
    residual_type Jr(2);
    jAction.applyJacobian(r, Jr);
    std::cout << *Jr.data() << std::endl;

    if (iterCount_==1)
    {
      // on this call, we know that
      // - correction should be all zeros
      // - nonlinear state should be the initial cond, i.e. [2,3]

      trueX(0) = 0.;
      trueX(1) = 0.;
      // this is because of the initial condition on the state
      trueR(0) = 4.;
      trueR(1) = 6.;
      trueJr(0) = 2.*4. + 3.*6.;
      trueJr(1) = 2.*4.*2. + 3.*6.*2.;
      if( !x.data()->isApprox(trueX) ) checkStr_="FAILED";
      if( !r.data()->isApprox(trueR) ) checkStr_="FAILED";
      if( !Jr.data()->isApprox(trueJr) ) checkStr_="FAILED";
    }

    if (iterCount_==2)
    {
      // on this call, we know that
      // - correction should be [-1,-1] because of sign convention in pressio
      // - nonlinear state should be [1,2]
      trueX(0) = -1.;
      trueX(1) = -1.;
      trueR(0) = 2.;
      trueR(1) = 4.;
      trueJr(0) = 1.*2. + 2.*4.;
      trueJr(1) = 1.*2.*2. + 2.*4.*2.;
      if( !x.data()->isApprox(trueX) ) checkStr_="FAILED";
      if( !r.data()->isApprox(trueR) ) checkStr_="FAILED";
      if( !Jr.data()->isApprox(trueJr) ) checkStr_="FAILED";
    }

    x(0) = 1.0;
    x(1) = 1.0;
  }
};

int main()
{
  pressio::log::initialize(pressio::logto::fileAndTerminal, "log.txt");
  pressio::log::setVerbosity({pressio::log::level::trace, pressio::log::level::trace});
  std::string strOut = "PASSED";

  using namespace pressio;
  using namespace pressio::solvers;
  using problem_t  = ValidSystem;
  using state_t	   = problem_t::state_type;

  // problem object
  problem_t sys;

  state_t state(2);
  state(0) = 2.; state(1) = 3.;

  // linear solver
  mylinsolver linearSolverObj(strOut);
  // nonlinear solver
  using pressio::solvers::nonlinear::experimental::createJacobianFreeNewtonRaphson;
  auto NonLinSolver = createJacobianFreeNewtonRaphson(sys, state, linearSolverObj);
  NonLinSolver.setMaxIterations(2);
  NonLinSolver.setStoppingCriterion(pressio::solvers::nonlinear::stop::afterMaxIters);
  NonLinSolver.solve(sys, state);

  const auto e1 = std::abs(state(0) - (1.));
  const auto e2 = std::abs(state(1) - (2.));
  if (e1>1e-13 or e2>1e-13) strOut = "FAILED";
  std::cout <<  strOut << std::endl;
  std::cout << *state.data() << std::endl;
}
