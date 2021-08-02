
#include "pressio_solvers.hpp"

struct CustomUpdate
{
  void reset(){}

  template<typename system_t, typename state_t, typename solver_mixin_t>
  void operator()(const system_t & sys,
		  state_t & state,
		  solver_mixin_t & solver)
  {
    PRESSIOLOG_DEBUG("nonlinsolver: custom update");
    const auto & correction = solver.correctionCRef();
    state(0) = correction(0) + 2.;
    state(1) = correction(1) + 3.;
  }
};

struct FakeProblem
{
  using scalar_type	= double;
  using state_type	= Eigen::VectorXd;
  using residual_type	= state_type;
  using jacobian_type	= Eigen::MatrixXd;

  residual_type createResidual() const{return residual_type(10);}
  jacobian_type createJacobian() const{return jacobian_type(10,2);}

  void residual(const state_type& x, residual_type & res) const
  {
    res.setConstant(1.);
  }

  void jacobian(const state_type & x, jacobian_type & jac) const
  {
    jac.setConstant(1.);
  }
};

template<class T>
struct FakeLinS
{
  int count_=0;
  using matrix_type = T;

  template<typename A_t, typename b_t, typename x_t>
  void solve(const A_t & A, const b_t & b, x_t & x)
  {
    ++count_;
    std::cout << x << std::endl;
    x.setConstant(1.1);
  }
};

int main()
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::trace});

  std::string sentinel = "PASSED";
  using namespace pressio;
  using problem_t = FakeProblem;
  using state_t	  = typename problem_t::state_type;
  using mat_t     = typename problem_t::jacobian_type;
  using hessian_t = mat_t;

  problem_t problem;
  state_t x(2);
  x(0)=1.0; x(1)=2.;

  FakeLinS<hessian_t> linSolver;

  auto GNSolver = nonlinearsolvers::createGaussNewton(problem,x,linSolver);
  auto criterion = nonlinearsolvers::stop::afterMaxIters;
  GNSolver.setStoppingCriterion(criterion);
  GNSolver.setMaxIterations(2);
  CustomUpdate U;
  GNSolver.solve(problem, x, U);

  if (std::abs(x(0)-3.1) > 1e-12){
    sentinel = "FAILED";
  }
  if (std::abs(x(1)-4.1) > 1e-12){
    sentinel = "FAILED";
  }
  std::cout << sentinel << std::endl;
  std::cout << std::setprecision(14) << x << std::endl;
  pressio::log::finalize();
  return 0;
}
