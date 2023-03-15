
#include <gtest/gtest.h>
#include "pressio/solvers_nonlinear.hpp"

using MyStateType = Eigen::VectorXd;
using MyGradientType = Eigen::VectorXd;
using MyJacType   = Eigen::SparseMatrix<double>;
using MyHessType  = Eigen::MatrixXd;

struct ValidSystemResJac
{
  using state_type     = MyStateType;
  using residual_type  = state_type;
  using jacobian_type  = MyJacType;

  state_type createState() const;
  residual_type createResidual() const;
  jacobian_type createJacobian() const;

  void residual(const state_type& x, residual_type & res) const;
  void jacobian(const state_type& x, jacobian_type & jac) const;
};

struct NonValidSystemResJac
{
  using state_type     = MyStateType;
  using residual_type  = state_type;
  using jacobian_type  = MyJacType;

  //residual_type createResidual() const; // missing on purpose to fail assert
  state_type createState() const;
  jacobian_type createJacobian() const;
  void residual(const state_type& x, residual_type & res) const;
  void jacobian(const state_type& x, jacobian_type & jac) const;
};

struct ValidSystemFusedResJac
{
  using state_type     = MyStateType;
  using residual_type  = state_type;
  using jacobian_type  = MyJacType;

  state_type createState() const;
  residual_type createResidual() const;
  jacobian_type createJacobian() const;

  void residualAndJacobian(const state_type& x,
         residual_type & res,
         jacobian_type & J,
         bool recomputeJacobian) const;
};

struct ValidSystemHessGrad {
  using state_type  = MyStateType;
  using hessian_type  = MyHessType;
  using gradient_type = MyGradientType;
  using residual_norm_type = double;

  state_type createState() const;
  hessian_type createHessian() const;
  gradient_type createGradient() const;

  void residualNorm(const state_type & state,
		    pressio::Norm normKind,
		    residual_norm_type & resNorm) const;

  void gradient(const state_type &,
		gradient_type &,
		::pressio::Norm normKind,
		double & normResidual,
		bool recomputeJacobian) const;

  void hessian(const state_type &, hessian_type &) const;
};

struct ValidSystemFusedHessGrad {
  using state_type  = MyStateType;
  using hessian_type  = MyHessType;
  using gradient_type = MyGradientType;
  using residual_norm_type = double;

  state_type createState() const;
  hessian_type createHessian() const;
  gradient_type createGradient() const;

  void residualNorm(const state_type & state,
		    pressio::Norm normKind,
		    residual_norm_type & resNorm) const;

  void hessianAndGradient(const state_type &, hessian_type &,
			  gradient_type &,
			  ::pressio::Norm normKind,
			  residual_norm_type & normResidual,
			  bool recomputeJacobian) const;
};

struct ValidLinearSolverNewtonRaphson
{
  using matrix_type = MyJacType;
  void solve(const matrix_type &, const MyStateType &, MyStateType &);
};

struct ValidLinearSolverGnOrLM
{
  using matrix_type = Eigen::MatrixXd;
  void solve(const matrix_type &, const MyStateType &, MyStateType &);
};

struct InvalidLinearSolver {
  using matrix_type = std::vector<std::vector<double>>;
};

template <typename A_t, typename r_t, typename state_type>
struct ValidQRSolver
{
  void computeThin(const A_t &);
  void applyQTranspose(const r_t &, state_type &) const;
  void applyRTranspose(const state_type &, state_type &) const;
  void solve(const state_type &, state_type &) const;
};

TEST(solvers_nonlinear, res_jac_api)
{
  using namespace pressio::nonlinearsolvers;
  using system_t   = ValidSystemResJac;
#ifdef PRESSIO_ENABLE_CXX20
  static_assert(SystemWithResidualAndJacobian<system_t>, "");
  static_assert(!SystemWithFusedResidualAndJacobian<system_t>, "");
#else
  static_assert(SystemWithResidualAndJacobian<system_t>::value, "");
  static_assert(!SystemWithFusedResidualAndJacobian<system_t>::value, "");
#endif
}

TEST(solvers_nonlinear, system_non_admissible_res_jac_api){
  using namespace pressio::nonlinearsolvers;
  using system_t   = NonValidSystemResJac;
#ifdef PRESSIO_ENABLE_CXX20
  static_assert(!SystemWithResidualAndJacobian<system_t>, "");
  static_assert(!SystemWithFusedResidualAndJacobian<system_t>, "");
#else
  static_assert(!SystemWithResidualAndJacobian<system_t>::value, "");
  static_assert(!SystemWithFusedResidualAndJacobian<system_t>::value, "");
#endif
}

TEST(solvers_nonlinear, system_admissible_fused_res_jac_api){
  using namespace pressio::nonlinearsolvers;
  using system_t   = ValidSystemFusedResJac;
#ifdef PRESSIO_ENABLE_CXX20
  static_assert(!SystemWithResidualAndJacobian<system_t>, "");
  static_assert(SystemWithFusedResidualAndJacobian<system_t>, "");
#else
  static_assert(!SystemWithResidualAndJacobian<system_t>::value, "");
  static_assert(SystemWithFusedResidualAndJacobian<system_t>::value, "");
#endif
}

// TEST(solvers_nonlinear, hes_gra_api){
//   using namespace pressio::nonlinearsolvers;
//   using system_t   = ValidSystemHessGrad;
// #ifdef PRESSIO_ENABLE_CXX20
//   static_assert(SystemWithHessianAndGradient<system_t>, "");
//   static_assert(!SystemWithFusedHessianAndGradient<system_t>, "");
// #else
//   static_assert(SystemWithHessianAndGradient<system_t>::value, "");
//   static_assert(!SystemWithFusedHessianAndGradient<system_t>::value, "");
// #endif
// }

// TEST(solvers_nonlinear, fused_hes_gra_api){
//   using namespace pressio::nonlinearsolvers;
//   using system_t   = ValidSystemFusedHessGrad;
// #ifdef PRESSIO_ENABLE_CXX20
//   static_assert(!SystemWithHessianAndGradient<system_t>, "");
//   static_assert(SystemWithFusedHessianAndGradient<system_t>, "");
// #else
//   static_assert(!SystemWithHessianAndGradient<system_t>::value, "");
//   static_assert(SystemWithFusedHessianAndGradient<system_t>::value, "");
// #endif
// }

TEST(solvers_meta, admissible_linear_solver_newtonraphon)
{
  using namespace pressio::nonlinearsolvers;
#ifdef PRESSIO_ENABLE_CXX20
  static_assert(LinearSolverForNewtonRaphson<
		ValidLinearSolverNewtonRaphson, MyJacType, MyStateType, MyStateType>, "");
  static_assert(!LinearSolverForNewtonRaphson<
		InvalidLinearSolver, MyJacType, MyStateType, MyStateType>, "");
#else
  static_assert(LinearSolverForNewtonRaphson<
		ValidLinearSolverNewtonRaphson, MyJacType, MyStateType, MyStateType>::value, "");
  static_assert(!LinearSolverForNewtonRaphson<
		InvalidLinearSolver, MyJacType, MyStateType, MyStateType>::value, "");
#endif
}

TEST(solvers_meta, admissible_linear_solver_nonlinear_ls)
{
  using state_type    = MyStateType;
  using namespace pressio::nonlinearsolvers;
#ifdef PRESSIO_ENABLE_CXX20
  static_assert(LinearSolverForNonlinearLeastSquares<ValidLinearSolverGnOrLM, state_type>, "");
  static_assert(!LinearSolverForNonlinearLeastSquares<InvalidLinearSolver, state_type>, "");
#else
  static_assert(LinearSolverForNonlinearLeastSquares<ValidLinearSolverGnOrLM, state_type>::value, "");
  static_assert(!LinearSolverForNonlinearLeastSquares<InvalidLinearSolver, state_type>::value, "");
#endif
}

TEST(solvers_nonlinear, admissible_qr_solver)
{
  using A_t = std::vector<MyStateType>;
  using r_t = MyStateType;
  using state_type = MyStateType;
  using solver_t = ValidQRSolver<A_t, r_t, state_type>;

#ifdef PRESSIO_ENABLE_CXX20
  static_assert(pressio::nonlinearsolvers::QRSolverForGnQr<solver_t, state_type, A_t, r_t>, "");
#else
  static_assert(pressio::nonlinearsolvers::QRSolverForGnQr<solver_t, state_type, A_t, r_t>::value, "");
#endif
}
