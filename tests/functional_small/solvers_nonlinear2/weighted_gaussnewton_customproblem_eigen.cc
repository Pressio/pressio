
#include "pressio/solvers_linear.hpp"
#include "pressio/solvers_nonlinear.hpp"
#include "./problems/problem3.hpp"

template<class scalar_t = double>
struct MyProblem
{
  using scalar_type	= scalar_t;
  using state_type	= Eigen::VectorXd;
  using residual_type	= state_type;
  using jacobian_type	= Eigen::MatrixXd;

  const double times_[8] = {1.,2.,3.,4., 5.,6.,7.,8};
  const double y_[8] = {3.29, 4.27, 5.3, 7.1, 10.1, 9.8, 16.1, 20.2};


  state_type createState() const {
    state_type a(2);
    a.setZero();
    return a;
  }

  residual_type createResidual() const {
    residual_type a(8);
    a.setZero();
    return a;
  }

  jacobian_type createJacobian() const {
    jacobian_type a(8, 2);
    a.setZero();
    return a;
  }

  void residual(const state_type& x, residual_type & res) const
  {
    for (auto i = 0; i < res.size(); i++) {
      res(i) = x(0) * exp(x(1)*times_[i]) - y_[i];
    }
  }

  void jacobian(const state_type & x, jacobian_type & jac) const {
    for (int i = 0; i < 8; i++) {
      double expval = exp(x(1) * times_[i]);
      jac(i,0) = expval;
      jac(i,1) = x(0)*times_[i]*expval;
    }
  }
};

template<class scalar_t = double>
struct Weigher{

  template<class T>
  void operator()(const Eigen::MatrixBase<T> & operand,
		  Eigen::MatrixBase<T> & result) const
  {}
};

struct MyLinSolver{
  void solve(const Eigen::MatrixXd & H,
	     const Eigen::VectorXd & b,
	     Eigen::VectorXd & x)
  {}
};

int main()
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::trace});

  std::string sentinel = "PASSED";
  using namespace pressio;
  using problem_t   = solvers::test::Problem3<double>;
  using state_t	    = typename problem_t::state_type;
  using mat_t   = typename problem_t::jacobian_type;
  using hessian_t = mat_t;

  problem_t problem;
  state_t x(2); x(0) = 2.0; x(1) = 0.25;

  auto GNSolver = create_gauss_newton_solver(problem,
					     MyLinSolver(),
					     Weigher{});

  std::cout << std::setprecision(14) << x << std::endl;
  pressio::log::finalize();
  return 0;
}
