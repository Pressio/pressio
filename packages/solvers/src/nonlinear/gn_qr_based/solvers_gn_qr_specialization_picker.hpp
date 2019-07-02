
#ifndef SOLVERS_GN_QR_SPECIALIZATION_PICKER_HPP_
#define SOLVERS_GN_QR_SPECIALIZATION_PICKER_HPP_

#include "../../solvers_fwd.hpp"
#include "../../meta/solvers_is_legitimate_system_for_nonlinear_solver.hpp"
#include "../../meta/solvers_is_legitimate_qr_solver_for_gn_qr.hpp"
#include "../../meta/solvers_is_non_default_line_search_tag.hpp"
#include "../../meta/solvers_is_non_default_convergence_tag.hpp"

namespace pressio{ namespace solvers{ namespace iterative{ namespace impl{

template <typename ... Args>
struct GNQRSpecializationPicker{

  // verify the sequence contains a valid system type
  using ic1 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::solvers::meta::is_legitimate_system_for_nonlinear_solver, Args...>;
  // store the type
  using system_t = ::pressio::mpl::variadic::at_or_t<void, ic1::value, Args...>;
  static_assert(!std::is_void<system_t>::value and
		ic1::value < sizeof... (Args),
		"A valid system type must be passed to GN templates");


  // verify the sequence contains a valid QR solver type
  using ic2 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::solvers::meta::is_legitimate_qr_solver_for_gn_qr, Args...>;
  // store the type
  using qr_solver_t = ::pressio::mpl::variadic::at_or_t<void, ic2::value, Args...>;
  static_assert(!std::is_void<qr_solver_t>::value and
		ic2::value < sizeof... (Args),
  		"A valid QR solver type must be passed to GN templates");
  // get the type of the matrix from the QR solver
  using qr_solver_matrix_t = typename ::pressio::qr::details::traits<qr_solver_t>::matrix_t;


  // check if sequence contains a non default line search tag
  using ic4 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::solvers::meta::is_non_default_line_search_tag, Args...>;
  // store the type
  using default_no_ls = ::pressio::solvers::iterative::gn::noLineSearch;
  using line_search_t = ::pressio::mpl::variadic::at_or_t<default_no_ls, ic4::value, Args...>;
  static_assert(!std::is_void<line_search_t>::value,
  		"The line search type for GN cannot be void");


  // check if sequence contains a non default convergence
  using ic5 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::solvers::meta::is_non_default_convergence_tag, Args...>;
  // store the type
  using default_conv = ::pressio::solvers::iterative::default_convergence;
  using convergence_t = ::pressio::mpl::variadic::at_or_t<default_conv, ic5::value, Args...>;
  static_assert(!std::is_void<convergence_t>::value,
  		"The convergence type for GN cannot be void");

  // using types detected above, define type of GN solver implementation
  using scalar_t = typename system_t::scalar_type;

  using type = ::pressio::solvers::iterative::impl::GaussNewtonQR<
    system_t, qr_solver_t, scalar_t, line_search_t, convergence_t>;
};

}}}}//end namespace pressio::solvers::iterative::impl
#endif
