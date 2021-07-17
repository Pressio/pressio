
#include <gtest/gtest.h>
#include "pressio_solvers.hpp"

struct ValidSolver {
  using matrix_type = Eigen::MatrixXd;

  template  <typename state_type>
  void solve(const matrix_type &, const state_type &, state_type &);
};

struct InvalidSolver {
  using matrix_type = Eigen::MatrixXd;
  // template  <typename state_type>
  // void solve(const matrix_type &, const state_type &, state_type &) const;
};

TEST(solvers_meta, admissible_linear_solver_newtonraphon)
{
  using state_type    = Eigen::VectorXd;
  using namespace pressio;
  static_assert(nonlinearsolvers::constraints::linear_solver_for_newton_raphson<ValidSolver, state_type>::value, "");
  static_assert(!nonlinearsolvers::constraints::linear_solver_for_newton_raphson<InvalidSolver, state_type>::value, "");
}

TEST(solvers_meta, admissible_linear_solver_nonlinear_ls)
{
  using state_type    = Eigen::VectorXd;
  using namespace pressio;
  static_assert(nonlinearsolvers::constraints::linear_solver_for_nonlinear_least_squares<ValidSolver, state_type>::value, "");
  static_assert(!nonlinearsolvers::constraints::linear_solver_for_nonlinear_least_squares<InvalidSolver, state_type>::value, "");
}
