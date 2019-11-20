
#include <gtest/gtest.h>
#include "CONTAINERS_ALL"
#include "SOLVERS_NONLINEAR"

struct System {
  using matrix_n_t = Eigen::SparseMatrix<double>;
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
	  typename T1, typename T2, typename T3,
	  typename T4=::pressio::solvers::iterative::gn::noLineSearch,
	  typename T5=::pressio::solvers::iterative::default_convergence,
	  typename T6=void>
struct checkTypes{

  using scalar_t		= typename pick_t::scalar_t;
  using system_t		= typename pick_t::system_t;
  using linear_solver_t		= typename pick_t::linear_solver_t;
  // using hessian_t		= typename pick_t::hessian_t;
  using line_search_t		= typename pick_t::line_search_t;
  using convergence_t		= typename pick_t::convergence_t;
  // using observer_when_conv_t	= typename pick_t::observer_when_conv_t;
  // using observer_each_step_t	= typename pick_t::observer_each_step_t;

  static_assert(std::is_same<scalar_t,		double>::value, "dont match");
  static_assert(std::is_same<system_t,		T1>::value, "dont match");
  static_assert(std::is_same<linear_solver_t,	T2>::value, "dont match");
  // static_assert(std::is_same<hessian_t,		T3>::value, "dont match");
  static_assert(std::is_same<line_search_t,	T4>::value, "dont match");
  static_assert(std::is_same<convergence_t,	T5>::value, "dont match");
  // static_assert(std::is_same<observer_when_conv_t, T6>::value, "dont match");
  // static_assert(std::is_same<observer_each_step_t, T6>::value, "dont match");
  static constexpr bool value = true;
};


TEST(solvers_meta, gn_normal_equations){
  using namespace pressio;
  using sys_t   = System;

  static_assert
    (solvers::meta::is_legitimate_system_for_gauss_newton_normal_eq
     <sys_t>::value, "");

  using hessian_type = containers::Matrix<Eigen::MatrixXd>;
  using tag = solvers::linear::iterative::LSCG;
  using lin_solver_t = solvers::iterative::EigenIterative<tag, hessian_type>;

  // define types, then rotate, it should not matter
  {
    using picker_t = solvers::iterative::impl::GNNEQSpecializationPicker<sys_t, lin_solver_t, hessian_type>;
    using gn_solver_t = solvers::iterative::GaussNewton<sys_t, lin_solver_t, hessian_type>;
    static_assert( checkTypes<picker_t, sys_t, lin_solver_t, hessian_type>::value, "");
    static_assert( !std::is_void<gn_solver_t>::value, "");
  }
  {
    using picker_t = solvers::iterative::impl::GNNEQSpecializationPicker<sys_t, hessian_type, lin_solver_t>;
    using gn_solver_t = solvers::iterative::GaussNewton<sys_t, hessian_type, lin_solver_t>;
    static_assert( checkTypes<picker_t, sys_t, lin_solver_t, hessian_type>::value, "");
    static_assert( !std::is_void<gn_solver_t>::value, "");
  }
  {
    using picker_t = solvers::iterative::impl::GNNEQSpecializationPicker<hessian_type, sys_t, lin_solver_t>;
    using gn_solver_t = solvers::iterative::GaussNewton<hessian_type, sys_t, lin_solver_t>;
    static_assert( checkTypes<picker_t, sys_t, lin_solver_t, hessian_type>::value, "");
    static_assert( !std::is_void<gn_solver_t>::value, "");
  }
  {
    using picker_t = solvers::iterative::impl::GNNEQSpecializationPicker<hessian_type, lin_solver_t, sys_t>;
    using gn_solver_t = solvers::iterative::GaussNewton<hessian_type, lin_solver_t, sys_t>;
    static_assert( checkTypes<picker_t, sys_t, lin_solver_t, hessian_type>::value, "");
    static_assert( !std::is_void<gn_solver_t>::value, "");
  }
}


TEST(solvers_meta, gn_normal_equations_noHess){
  using namespace pressio;
  using sys_t   = System;

  static_assert
    (solvers::meta::is_legitimate_system_for_gauss_newton_normal_eq
     <sys_t>::value, "");

  using hessian_type = containers::Matrix<Eigen::MatrixXd>;
  using tag = solvers::linear::iterative::LSCG;
  using lin_solver_t = solvers::iterative::EigenIterative<tag, hessian_type>;

  // define types, then rotate, it should not matter
  {
    using picker_t = solvers::iterative::impl::GNNEQSpecializationPicker<sys_t, lin_solver_t>;
    using gn_solver_t = solvers::iterative::GaussNewton<sys_t, lin_solver_t>;
    static_assert( checkTypes<picker_t, sys_t, lin_solver_t, hessian_type>::value, "");
    static_assert( !std::is_void<gn_solver_t>::value, "");
  }
  {
    using picker_t = solvers::iterative::impl::GNNEQSpecializationPicker<sys_t, lin_solver_t>;
    using gn_solver_t = solvers::iterative::GaussNewton<sys_t, lin_solver_t>;
    static_assert( checkTypes<picker_t, sys_t, lin_solver_t, hessian_type>::value, "");
    static_assert( !std::is_void<gn_solver_t>::value, "");
  }
  {
    using picker_t = solvers::iterative::impl::GNNEQSpecializationPicker<sys_t, lin_solver_t>;
    using gn_solver_t = solvers::iterative::GaussNewton<sys_t, lin_solver_t>;
    static_assert( checkTypes<picker_t, sys_t, lin_solver_t, hessian_type>::value, "");
    static_assert( !std::is_void<gn_solver_t>::value, "");
  }
  {
    using picker_t = solvers::iterative::impl::GNNEQSpecializationPicker<lin_solver_t, sys_t>;
    using gn_solver_t = solvers::iterative::GaussNewton<lin_solver_t, sys_t>;
    static_assert( checkTypes<picker_t, sys_t, lin_solver_t, hessian_type>::value, "");
    static_assert( !std::is_void<gn_solver_t>::value, "");
  }
}


TEST(solvers_meta, gn_normal_equations_nondef_conv){
  using namespace pressio;
  using sys_t   = System;

  static_assert
    (solvers::meta::is_legitimate_system_for_gauss_newton_normal_eq
     <sys_t>::value, "");

  using hessian_type = containers::Matrix<Eigen::MatrixXd>;
  using tag = solvers::linear::iterative::LSCG;
  using lin_solver_t = solvers::iterative::EigenIterative<tag, hessian_type>;

  using ls_t   = solvers::iterative::gn::ArmijoLineSearch;
  using conv_t = solvers::iterative::converged_when::absoluteNormCorrectionBelowTol;

  // define types, then rotate, it should not matter
  {
    using picker_t = solvers::iterative::impl::GNNEQSpecializationPicker<sys_t, lin_solver_t, hessian_type, ls_t, conv_t>;
    using gn_solver_t = solvers::iterative::GaussNewton<sys_t, lin_solver_t, hessian_type, ls_t, conv_t>;
    static_assert( checkTypes<picker_t, sys_t, lin_solver_t, hessian_type, ls_t, conv_t>::value, "");
    static_assert( !std::is_void<gn_solver_t>::value, "");
  }
  {
    using picker_t = solvers::iterative::impl::GNNEQSpecializationPicker<sys_t, hessian_type, ls_t, conv_t, lin_solver_t>;
    using gn_solver_t = solvers::iterative::GaussNewton<sys_t, hessian_type, ls_t, conv_t, lin_solver_t>;
    static_assert( checkTypes<picker_t, sys_t, lin_solver_t, hessian_type, ls_t, conv_t>::value, "");
    static_assert( !std::is_void<gn_solver_t>::value, "");
  }
  {
    using picker_t = solvers::iterative::impl::GNNEQSpecializationPicker<hessian_type, conv_t, ls_t, sys_t, lin_solver_t>;
    using gn_solver_t = solvers::iterative::GaussNewton<hessian_type, conv_t, ls_t, sys_t, lin_solver_t>;
    static_assert( checkTypes<picker_t, sys_t, lin_solver_t, hessian_type, ls_t, conv_t>::value, "");
    static_assert( !std::is_void<gn_solver_t>::value, "");
  }
  {
    using picker_t = solvers::iterative::impl::GNNEQSpecializationPicker<ls_t, conv_t, hessian_type, lin_solver_t, sys_t>;
    using gn_solver_t = solvers::iterative::GaussNewton<ls_t, conv_t, hessian_type, lin_solver_t, sys_t>;
    static_assert( checkTypes<picker_t, sys_t, lin_solver_t, hessian_type, ls_t, conv_t>::value, "");
    static_assert( !std::is_void<gn_solver_t>::value, "");
  }
}
