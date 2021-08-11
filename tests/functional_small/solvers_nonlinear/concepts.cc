
#include <gtest/gtest.h>
#include "pressio/solvers_nonlinear.hpp"

struct ValidSystemResJac 
{
  using scalar_type    = double;
  using state_type     = std::vector<scalar_type>;
  using residual_type  = state_type;
  using jacobian_type  = std::vector<std::vector<scalar_type>>;

  residual_type createResidual() const;
  jacobian_type createJacobian() const;

  void residual(const state_type& x, residual_type & res) const;
  void jacobian(const state_type& x, jacobian_type & jac) const;
};

struct NonValidSystemResJac 
{
  using scalar_type     = double;
  using state_type     = std::vector<scalar_type>;
  using residual_type  = state_type;
  using jacobian_type  = std::vector<std::vector<scalar_type>>;

  //residual_type createResidual() const; // missing on purpose to fail assert
  jacobian_type createJacobian() const;
  void residual(const state_type& x, residual_type & res) const;
  void jacobian(const state_type& x, jacobian_type & jac) const;
};

struct ValidSystemFusedResJac
{
  using scalar_type     = double;
  using state_type     = std::vector<scalar_type>;
  using residual_type  = state_type;
  using jacobian_type  = std::vector<std::vector<scalar_type>>;

  residual_type createResidual() const;
  jacobian_type createJacobian() const;

  void residualAndJacobian(const state_type& x,
         residual_type & res, 
         jacobian_type & J,
         bool recomputeJacobian) const;
};

struct ValidSystemHessGrad {
  using scalar_type     = double;
  using state_type  = std::vector<scalar_type>;
  using hessian_type  = std::vector<std::vector<scalar_type>>;
  using gradient_type = state_type;

  hessian_type createHessian() const;
  gradient_type createGradient() const;

  void residualNorm(const state_type & state,
          pressio::Norm normKind,
          scalar_type & resNorm) const;

  void gradient(const state_type &,
    gradient_type &,
    ::pressio::Norm normKind,
    scalar_type & normResidual,
    bool recomputeJacobian) const;

  void hessian(const state_type &, hessian_type &) const;
};

struct ValidSystemFusedHessGrad {
  using scalar_type     = double;
  using state_type  = std::vector<scalar_type>;
  using hessian_type  = std::vector<std::vector<scalar_type>>;
  using gradient_type = state_type;

  hessian_type createHessian() const;
  gradient_type createGradient() const;

  void residualNorm(const state_type & state,
          pressio::Norm normKind,
          scalar_type & resNorm) const;

  void hessianAndGradient(const state_type &, hessian_type &,
        gradient_type &,
        ::pressio::Norm normKind,
        scalar_type & normResidual,
        bool recomputeJacobian) const;
};

TEST(solvers_nonlinear, res_jac_api)
{
  using namespace pressio;
  using system_t   = ValidSystemResJac;
  static_assert(nonlinearsolvers::compliant_with_residual_jacobian_api<system_t>::value, "");
  static_assert(!nonlinearsolvers::compliant_with_fused_residual_jacobian_api<system_t>::value, "");
}

TEST(solvers_nonlinear, system_non_admissible_res_jac_api){
  using namespace pressio;
  using system_t   = NonValidSystemResJac;
  static_assert(!nonlinearsolvers::compliant_with_residual_jacobian_api<system_t>::value, "");
  static_assert(!nonlinearsolvers::compliant_with_fused_residual_jacobian_api<system_t>::value, "");
}

TEST(solvers_nonlinear, system_admissible_fused_res_jac_api){
  using namespace pressio;
  using system_t   = ValidSystemFusedResJac;
  static_assert(!nonlinearsolvers::compliant_with_residual_jacobian_api<system_t>::value, "");
  static_assert(nonlinearsolvers::compliant_with_fused_residual_jacobian_api<system_t>::value, "");
}

TEST(solvers_nonlinear, hes_gra_api){
  using namespace pressio;
  using system_t   = ValidSystemHessGrad;
  static_assert(nonlinearsolvers::compliant_with_hessian_gradient_api<system_t>::value, "");
  static_assert(!nonlinearsolvers::compliant_with_fused_hessian_gradient_api<system_t>::value, "");
}


TEST(solvers_nonlinear, fused_hes_gra_api){
  using namespace pressio;
  using system_t   = ValidSystemFusedHessGrad;
  static_assert(!nonlinearsolvers::compliant_with_hessian_gradient_api<system_t>::value, "");
  static_assert(nonlinearsolvers::compliant_with_fused_hessian_gradient_api<system_t>::value, "");
}

struct ValidLinearSolver 
{
  using matrix_type = std::vector<std::vector<double>>;

  template  <typename state_type>
  void solve(const matrix_type &, const state_type &, state_type &);
};

struct InvalidLinearSolver {
  using matrix_type = std::vector<std::vector<double>>;

  // template  <typename state_type>
  // void solve(const matrix_type &, const state_type &, state_type &) const;
};


TEST(solvers_meta, admissible_linear_solver_newtonraphon)
{
  using state_type    = std::vector<double>;
  using namespace pressio;
  static_assert(nonlinearsolvers::admissible_linear_solver_for_newton_raphson<ValidLinearSolver, state_type>::value, "");
  static_assert(!nonlinearsolvers::admissible_linear_solver_for_newton_raphson<InvalidLinearSolver, state_type>::value, "");
}

TEST(solvers_meta, admissible_linear_solver_nonlinear_ls)
{
  using state_type    = std::vector<double>;
  using namespace pressio;
  static_assert(nonlinearsolvers::admissible_linear_solver_for_nonlinear_least_squares<ValidLinearSolver, state_type>::value, "");
  static_assert(!nonlinearsolvers::admissible_linear_solver_for_nonlinear_least_squares<InvalidLinearSolver, state_type>::value, "");
}

template <typename A_t, typename r_t, typename state_type>
struct ValidQRSolver 
{
  void computeThin(const A_t &);
  void applyQTranspose(const r_t &, state_type &) const;
  void applyRTranspose(const state_type &, state_type &) const;
  void solve(const state_type &, state_type &) const;
};

TEST(solvers_nonlinear, admissible_qr_solver)
{
  using A_t = std::vector<std::vector<double>>;
  using r_t = std::vector<double>;
  using state_type = std::vector<double>;

  using solver_t = ValidQRSolver<A_t, r_t, state_type>;
  static_assert(pressio::nonlinearsolvers::admissible_qr_solver_for_gn_qr<
    solver_t, state_type, A_t, r_t>::value, "");
}
