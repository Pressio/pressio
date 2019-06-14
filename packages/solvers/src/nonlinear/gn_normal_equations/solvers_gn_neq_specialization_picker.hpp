
#ifndef SOLVERS_GN_NEQ_SPECIALIZATION_PICKER_HPP_
#define SOLVERS_GN_NEQ_SPECIALIZATION_PICKER_HPP_

#include "../../solvers_forward_declarations.hpp"
#include "../../meta/solvers_is_legitimate_system_for_nonlinear_solver.hpp"
#include "../../meta/solvers_is_legitimate_linear_solver_for_gn_normeq.hpp"
#include "../../meta/solvers_is_legitimate_hessian_for_gn_normeq.hpp"
#include "../../meta/solvers_is_non_default_line_search_tag.hpp"
#include "../../meta/solvers_is_non_default_convergence_tag.hpp"
#include "../../meta/solvers_is_legitimate_residual_observer_when_converged.hpp"

namespace rompp{ namespace solvers{ namespace iterative{ namespace impl{

template <typename T1, typename T2>
struct ObserverTypesSupported{
  static constexpr bool value = false;
};

template <>
struct ObserverTypesSupported<void, void>{
  static constexpr bool value = true;
  using type = void;
};

template <typename T>
struct ObserverTypesSupported<void, T>{
  static constexpr bool value = true;
  using type = T;
};

template <typename T>
struct ObserverTypesSupported<T, void>{
  static constexpr bool value = true;
  using type = T;
};

template <typename T>
struct ObserverTypesSupported<T, T>{
  static constexpr bool value = true;
  using type = T;
};




template <typename ... Args>
struct GNNEQSpecializationPicker{

  // verify the sequence contains a valid system type
  using ic1 = ::rompp::mpl::variadic::find_if_unary_pred_t<
    ::rompp::solvers::meta::is_legitimate_system_for_nonlinear_solver, Args...>;
  // store the type
  using system_t = ::rompp::mpl::variadic::at_or_t<void, ic1::value, Args...>;
  static_assert(!std::is_void<system_t>::value and
		ic1::value < sizeof... (Args),
		"A valid system type must be passed to GN templates");


  // verify the sequence contains a valid solver type
  using ic2 = ::rompp::mpl::variadic::find_if_unary_pred_t<
    ::rompp::solvers::meta::is_legitimate_linear_solver_for_gn_normeq, Args...>;
  // store the type
  using linear_solver_t = ::rompp::mpl::variadic::at_or_t<void, ic2::value, Args...>;
  static_assert(!std::is_void<linear_solver_t>::value and
		ic2::value < sizeof... (Args),
  		"A valid linear solver type must be passed to GN templates");
  // store the type of the matrix in the linear solver
  using linear_solver_matrix_t = typename linear_solver_t::matrix_type;


  // check if the sequence contains a valid hessian type
  using ic3 = ::rompp::mpl::variadic::find_if_unary_pred_t<
    ::rompp::solvers::meta::is_legitimate_hessian_for_gn_normeq, Args...>;
  // if no hessian is found, then get it from the linear solver type
  using hessian_t = ::rompp::mpl::variadic::at_or_t<linear_solver_matrix_t, ic3::value, Args...>;
  // in every scenario, the hesssian type must match the matrix type of linear solver
  static_assert( std::is_same<hessian_t, linear_solver_matrix_t>::value,
		 "Hessian type passed to GN must match the matrix type in the linear solver");
  static_assert(!std::is_void<hessian_t>::value,
  		"The hessian type cannot be void");


  // check if sequence contains a non default line search tag
  using ic4 = ::rompp::mpl::variadic::find_if_unary_pred_t<
    ::rompp::solvers::meta::is_non_default_line_search_tag, Args...>;
  // store the type
  using default_no_ls = ::rompp::solvers::iterative::gn::noLineSearch;
  using line_search_t = ::rompp::mpl::variadic::at_or_t<default_no_ls, ic4::value, Args...>;
  static_assert(!std::is_void<line_search_t>::value,
		"The line search type for GN cannot be void");


  // check if sequence contains a non default convergence
  using ic5 = ::rompp::mpl::variadic::find_if_unary_pred_t<
    ::rompp::solvers::meta::is_non_default_convergence_tag, Args...>;
  // store the type
  using default_conv = ::rompp::solvers::iterative::default_convergence;
  using convergence_t = ::rompp::mpl::variadic::at_or_t<default_conv, ic5::value, Args...>;
  static_assert(!std::is_void<convergence_t>::value,
		"The convergence type for GN cannot be void");


  // check if sequence contains an observer for the residual after GN is converged
  using ic6 = ::rompp::mpl::variadic::find_if_binary_pred_t<
    typename system_t::residual_type,
    ::rompp::solvers::meta::is_legitimate_residual_observer_when_solver_converged,
    Args...>;
  // store the type
  using observer_when_conv_t = ::rompp::mpl::variadic::at_or_t<void, ic6::value, Args...>;


  // check if sequence contains an observer for the residual at each GN step
  using ic7 = ::rompp::mpl::variadic::find_if_binary_pred_t<
    typename system_t::residual_type,
    ::rompp::solvers::meta::is_legitimate_residual_observer_each_solver_step, Args...>;
  // store the type
  using observer_each_step_t = ::rompp::mpl::variadic::at_or_t<void, ic7::value, Args...>;

  /* currently we only allows all observation methods on the residual to be
   * in the same observer class. So we need to check here that when both
   *  observers type detected are NOT void, then they must be the same type */
  using obs_supported_t = ObserverTypesSupported<observer_when_conv_t, observer_each_step_t>;
  static_assert( obs_supported_t::value,
  "Currently, can only accept a single residual observer type to do everything. \
  You can still have multiple methods inside that type to cover all possible observations,\
  but cannot pass multiple types with individual observation methods");
  using observer_t = typename obs_supported_t::type;


  // using types detected above, define type of GN solver implementation
  using scalar_t = typename system_t::scalar_type;

  using type = ::rompp::solvers::iterative::impl::GaussNewton<
    system_t, hessian_t, linear_solver_t, scalar_t, line_search_t,
    convergence_t, observer_t>;
};

}}}}//end namespace rompp::solvers::iterative::impl
#endif
