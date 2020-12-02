
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
    res(0) =  x(0)*x(0)*x(0) + x(1) - 1.0;
    res(1) = -x(0) + x(1)*x(1)*x(1) + 1.0;
  }

  void applyJacobian(const residual_type & in,
		     residual_type & out) const
  {
    jacobian_type jac(2,2);
    //computeJacobian(x, jac);
    *out.data() = (*jac.data()) * (*in.data());
  }

private:
  void computeJacobian(const state_type & x,
		       jacobian_type & jac) const
  {
    jac.data()->coeffRef(0, 0) = 3.0*x(0)*x(0);
    jac.data()->coeffRef(0, 1) =  1.0;
    jac.data()->coeffRef(1, 0) = -1.0;
    jac.data()->coeffRef(1, 1) = 3.0*x(1)*x(1);
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

  // solve Ax = b
  // where A only evaluates action of the jacobian
  void solve(const ValidSystem & A,
             const residual_type & b,
             state_type & x)
  {
    ++iterCount_;
    std::cout << iterCount_ << std::endl;
    std::cout << *x.data() << std::endl;
    std::cout << *b.data() << std::endl;

    eig_v_t trueB(2);
    eig_v_t trueX(2);
    eig_v_t trueJb(2);
    // residual_type Jb(2);
    // F.applyJacobian(x, b, Jb);
    // std::cout << *Jb.data() << std::endl;

    if (iterCount_==1){
      trueX(0) = 0.;
      trueX(1) = 0.;
      trueB(0) = 0.001*0.001*0.001 + 0.0001 - 1.;
      trueB(1) = -0.001 + 0.0001*0.0001*0.0001 + 1.;
      // trueJb(0) = 3.*0.001*0.001*trueB(0)+1.*trueB(1);
      // trueJb(1) = -1.*trueB(0)+1.*trueB(1);
      //if( !x.data()->isApprox(trueX) ) checkStr_="FAILED";
      // if( !b.data()->isApprox(trueB) ) checkStr_="FAILED";
      // if( !Jb.data()->isApprox(trueJb) ) checkStr_="FAILED";
    }

    // if (iterCount_==2){
    //   trueX(0) = -1.;
    //   trueX(1) = -1.;
    //   trueB(0) = 0.;
    //   trueB(1) = 0.9993;
    //   // if( !x.data()->isApprox(trueX) ) checkStr_="FAILED";
    //   // if( !b.data()->isApprox(trueB) ) checkStr_="FAILED";
    // }

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

  problem_t sys;
  state_t state(2);
  state(0) = 0.001; state(1) = 0.0001;

  // linear solver
  mylinsolver linearSolverObj(strOut);
  // nonlinear solver
  using pressio::solvers::nonlinear::experimental::createJacobianFreeNewtonRaphson;
  auto NonLinSolver = createJacobianFreeNewtonRaphson(sys, state, linearSolverObj);
  NonLinSolver.setMaxIterations(2);
  NonLinSolver.setStoppingCriterion(pressio::solvers::nonlinear::stop::afterMaxIters);
  NonLinSolver.solve(sys, state);

  // const auto e1 = std::abs(state(0) - (1.));
  // const auto e2 = std::abs(state(1) - (0.));
  // if (e1>1e-8 or e2>1e-8) strOut = "FAILED"
  std::cout <<  strOut << std::endl;
  std::cout << *state.data() << std::endl;
}
