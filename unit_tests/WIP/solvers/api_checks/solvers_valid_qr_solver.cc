
#include <gtest/gtest.h>
#include "pressio_solvers.hpp"

template <typename A_t, typename r_t, typename state_type>
struct ValidSolver 
{
  void computeThin(const A_t &);
  void applyQTranspose(const r_t &, state_type &) const;
  void applyRTranspose(const state_type &, state_type &) const;
  void solve(const state_type &, state_type &) const;
};

TEST(solvers_meta, admissible_qr_solver)
{
  using A_t = Eigen::MatrixXd;
  using r_t = Eigen::VectorXd;
  using state_type = Eigen::VectorXd;

  using solver_t = ValidSolver<A_t, r_t, state_type>;
  static_assert(pressio::nonlinearsolvers::constraints::qr_solver_for_gn_qr<
    solver_t, state_type, A_t, r_t>::value, "");

  using qr_solver_type = pressio::qr::QRSolver<A_t, pressio::qr::Householder>;
  static_assert(pressio::nonlinearsolvers::constraints::qr_solver_for_gn_qr<
    qr_solver_type, state_type, A_t, r_t>::value, "");
}
