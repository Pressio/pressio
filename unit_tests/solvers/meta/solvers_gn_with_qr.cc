
#include <gtest/gtest.h>
#include "CONTAINERS_ALL"
#include "SOLVERS_NONLINEAR"

struct System {
  using matrix_n_t = Eigen::MatrixXd;
  using matrix_w_t = pressio::containers::Matrix<matrix_n_t>;
  using vector_n_t = Eigen::VectorXd;
  using vector_w_t = pressio::containers::Vector<vector_n_t>;

  using scalar_type     = double;
  using state_type	= vector_w_t;
  using residual_type	= state_type;
  using jacobian_type	= matrix_w_t;

  void residual(const state_type & x, residual_type & res);
  residual_type residual(const state_type & x);

  void jacobian(const state_type & x, jacobian_type & jac);
  jacobian_type jacobian(const state_type & x);
};

template <typename pick_t,
	  typename T1, typename T2,
	  typename T3=::pressio::solvers::iterative::gn::noLineSearch,
	  typename T4=::pressio::solvers::iterative::default_convergence>
struct checkTypes{

  using scalar_t		= typename pick_t::scalar_t;
  using system_t		= typename pick_t::system_t;
  using solver_t		= typename pick_t::qr_solver_t;
  using line_search_t		= typename pick_t::line_search_t;
  using convergence_t		= typename pick_t::convergence_t;

  static_assert(std::is_same<scalar_t,		double>::value, "dont match");
  static_assert(std::is_same<system_t,		T1>::value, "dont match");
  static_assert(std::is_same<solver_t,		T2>::value, "dont match");
  static_assert(std::is_same<line_search_t,	T3>::value, "dont match");
  static_assert(std::is_same<convergence_t,	T4>::value, "dont match");
  static constexpr bool value = true;
};

TEST(solvers_meta, gn_with_qr){
  using namespace pressio;
  using sys_t   = System;

  static_assert
    (solvers::meta::is_legitimate_system_for_nonlinear_solver
     <sys_t>::value, "");

  using mat_type    = typename sys_t::jacobian_type;
  using qr_algo	    = qr::Householder;
  using qr_solver_t = qr::QRSolver<mat_type, qr_algo>;

  // define types, then rotate, it should not matter
  {
    static_assert(::pressio::solvers::meta::is_legitimate_qr_solver_for_gn_qr<qr_solver_t>::value, "" );

    using picker_t = solvers::iterative::impl::GNQRSpecializationPicker<sys_t, qr_solver_t>;
    using gn_solver_t = solvers::iterative::GaussNewtonQR<sys_t, qr_solver_t>;
    static_assert( checkTypes<picker_t, sys_t, qr_solver_t>::value, "");
    static_assert( !std::is_void<gn_solver_t>::value, "");
  }
  {
    using picker_t = solvers::iterative::impl::GNQRSpecializationPicker<qr_solver_t, sys_t>;
    using gn_solver_t = solvers::iterative::GaussNewtonQR<qr_solver_t, sys_t>;
    static_assert( checkTypes<picker_t, sys_t, qr_solver_t>::value, "");
    static_assert( !std::is_void<gn_solver_t>::value, "");
  }
}


TEST(solvers_meta, gn_normal_equations_nondef_conv){
  using namespace pressio;
  using sys_t   = System;

  static_assert
    (solvers::meta::is_legitimate_system_for_nonlinear_solver
     <sys_t>::value, "");

  using mat_type    = typename sys_t::jacobian_type;
  using qr_algo	    = qr::Householder;
  using qr_solver_t = qr::QRSolver<mat_type, qr_algo>;

  using ls_t   = solvers::iterative::gn::ArmijoLineSearch;
  using conv_t = solvers::iterative::converged_when::absoluteNormCorrectionBelowTol;

  // define types, then rotate, it should not matter
  {
    using picker_t = solvers::iterative::impl::GNQRSpecializationPicker<sys_t,qr_solver_t, ls_t, conv_t>;
    using gn_solver_t = solvers::iterative::GaussNewtonQR<sys_t,qr_solver_t,  ls_t, conv_t>;
    static_assert( checkTypes<picker_t, sys_t,qr_solver_t, ls_t, conv_t>::value, "");
    static_assert( !std::is_void<gn_solver_t>::value, "");
  }
  {
    using picker_t = solvers::iterative::impl::GNQRSpecializationPicker<sys_t,  ls_t, conv_t,qr_solver_t>;
    using gn_solver_t = solvers::iterative::GaussNewtonQR<sys_t,  ls_t, conv_t,qr_solver_t>;
    static_assert( checkTypes<picker_t, sys_t,qr_solver_t, ls_t, conv_t>::value, "");
    static_assert( !std::is_void<gn_solver_t>::value, "");
  }
  {
    using picker_t = solvers::iterative::impl::GNQRSpecializationPicker< conv_t, ls_t, sys_t,qr_solver_t>;
    using gn_solver_t = solvers::iterative::GaussNewtonQR< conv_t, ls_t, sys_t,qr_solver_t>;
    static_assert( checkTypes<picker_t, sys_t,qr_solver_t, ls_t, conv_t>::value, "");
    static_assert( !std::is_void<gn_solver_t>::value, "");
  }
  {
    using picker_t = solvers::iterative::impl::GNQRSpecializationPicker<ls_t, conv_t, qr_solver_t, sys_t>;
    using gn_solver_t = solvers::iterative::GaussNewtonQR<ls_t, conv_t, qr_solver_t, sys_t>;
    static_assert( checkTypes<picker_t, sys_t,qr_solver_t, ls_t, conv_t>::value, "");
    static_assert( !std::is_void<gn_solver_t>::value, "");
  }
}
