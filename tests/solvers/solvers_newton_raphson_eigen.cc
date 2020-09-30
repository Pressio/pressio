
#include "pressio_solvers.hpp"

struct ValidSystem {
  // Matrix typedefs
  using matrix_n_t = Eigen::SparseMatrix<double>;
  using matrix_w_t = pressio::containers::SparseMatrix<matrix_n_t>;

  // Vector typedefs
  using vector_n_t = Eigen::VectorXd;
  using vector_w_t = pressio::containers::Vector<vector_n_t>;

  using state_type = vector_w_t;
  using residual_type = state_type;
  using jacobian_type = matrix_w_t;
  using scalar_type = double;

  residual_type createResidual() const {
    return residual_type(2);
  }

  jacobian_type createJacobian() const {
    return jacobian_type(2, 2);
  }

  void residual(const state_type& x, 
                residual_type& res,
                pressio::Norm normKind, 
                scalar_type & norm) const 
  {
    res[0] =  x[0]*x[0]*x[0] + x[1] - 1.0;
    res[1] = -x[0] + x[1]*x[1]*x[1] + 1.0;
    if (normKind == pressio::Norm::L2) 
      norm = res.data()->norm();
    else if (normKind == pressio::Norm::L1) 
      norm = res.data()->lpNorm<1>();
  }

  void jacobian(const state_type& x, jacobian_type& jac) const {
    jac.data()->coeffRef(0, 0) = 3.0*x[0]*x[0];
    jac.data()->coeffRef(0, 1) =  1.0;
    jac.data()->coeffRef(1, 0) = -1.0;
    jac.data()->coeffRef(1, 1) = 3.0*x[1]*x[1];
  }
};


int main()
{
  using namespace pressio;
  using namespace pressio::solvers;

  using problem_t  = ValidSystem;
  using state_t	   = problem_t::state_type;
  using jacobian_t = problem_t::jacobian_type;

  problem_t sys;
  // my solution vector
  state_t y(2);
  // initial condition
  y[0] = 0.001; y[1] = 0.0001;

  // linear system
  using lin_solver_t = linear::Solver<linear::iterative::LSCG, jacobian_t>;
  lin_solver_t linearSolverObj;

  auto NonLinSolver = pressio::solvers::nonlinear::createNewtonRaphson(sys, y, linearSolverObj);

  NonLinSolver.solve(sys, y);

  std::string strOut = "PASSED";
  const auto e1 = std::abs(y[0] - (1.)); 
  const auto e2 = std::abs(y[1] - (0.)); 
  if (e1>1e-8 or e2>1e-8) strOut = "FAILED";

  std::cout <<  strOut << std::endl;
  std::cout << *y.data() << std::endl;
}
