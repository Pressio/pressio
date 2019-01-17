
#include <gtest/gtest.h>
#include "CORE_ALL"
#include "SOLVERS_NONLINEAR"
#include "QR_BASIC"

struct Rosenbrock {

  using eig_dyn_mat	= Eigen::MatrixXd;
  using eig_dyn_vec	= Eigen::VectorXd;
  using jacobian_w_t	= rompp::core::Matrix<eig_dyn_mat>;
  using state_w_t	= rompp::core::Vector<eig_dyn_vec>;
  using state_type	= state_w_t;
  using residual_type	= state_type;
  using jacobian_type	= jacobian_w_t;
  static constexpr int nf = 4; // num functions
  static constexpr int nv = 3; // num variables

  void residual(const state_type& x, residual_type & res) const {
    auto x1 = x[0];
    auto x2 = x[1];
    auto x3 = x[2];

    res[0] = 10.*(x3 - x2*x2);
    res[1] = 10.*(x2 - x1*x1);
    res[2] = (1.-x1);
    res[3] = (1.-x2);
  }

  residual_type residual(const state_type& x) const {
    residual_type res(nf);
    this->residual(x, res);
    return res;
  }

  void jacobian(const state_type & x, jacobian_type & jac) const {
    auto x1 = x[0];
    auto x2 = x[1];
    auto & JJ = *jac.data();
    JJ.setZero();

    JJ(0,1) = -20.*x2;
    JJ(0,2) = 10.;

    JJ(1,0) = -20.*x1;
    JJ(1,1) = 10.;

    JJ(2,0) = -1.;
    JJ(3,1) = -1.;
  }

  jacobian_type jacobian(const state_type& x) const {
    jacobian_type jac(nf, nv);
    this->jacobian(x, jac);
    return jac;
  }
};


TEST(solvers_nonlinear_least_squares, rosenbrockGNQRLineSearch){
  using namespace rompp;
  using state_w_t = core::Vector<Eigen::VectorXd>;
  using sc_t	  = double;
  using mat_type  = typename Rosenbrock::jacobian_w_t;
  Rosenbrock problem;

  state_w_t x(3);
  x[0] = -1.5; x[1] = 1.1; x[2] = 1.2;

  // define type of QR and GaussNewton solver
  using qr_algo = qr::Householder;
  using qr_type = qr::QRSolver<mat_type, qr_algo>;
  using lsearch_t = solvers::iterative::gn::ArmijoLineSearch;
  solvers::iterative::GaussNewtonQRLineSearch<sc_t,
					      qr_type,
					      lsearch_t> solver;

  solver.setTolerance(1e-8);
  solver.solve(problem, x);
  std::cout << std::setprecision(14) << *x.data() << std::endl;
  EXPECT_NEAR( x(0), 1., 1e-6 );
  EXPECT_NEAR( x(1), 1., 1e-6 );
  EXPECT_NEAR( x(2), 1., 1e-6 );
}

TEST(solvers_nonlinear_least_squares, rosenbrockGNQRLineSearchPassTypes){
  using namespace rompp;
  using state_w_t = core::Vector<Eigen::VectorXd>;
  using sc_t	  = double;
  using mat_type  = typename Rosenbrock::jacobian_w_t;
  using problem_t = Rosenbrock;
  problem_t problem;

  state_w_t x(3);
  x[0] = -1.5; x[1] = 1.1; x[2] = 1.2;

  // define type of QR and GaussNewton solver
  using qr_algo = qr::Householder;
  using qr_type = qr::QRSolver<mat_type, qr_algo>;
  using lsearch_t = solvers::iterative::gn::ArmijoLineSearch;
  using converged_when_t = solvers::iterative::default_convergence;
  using gnsolver_t = solvers::iterative::GaussNewtonQRLineSearch<
    sc_t, qr_type, lsearch_t, converged_when_t, problem_t>;
  gnsolver_t solver(problem, x);

  solver.setTolerance(1e-8);
  solver.solve(problem, x);

  std::cout << std::setprecision(14) << *x.data() << std::endl;
  EXPECT_NEAR( x(0), 1., 1e-6 );
  EXPECT_NEAR( x(1), 1., 1e-6 );
  EXPECT_NEAR( x(2), 1., 1e-6 );
}
