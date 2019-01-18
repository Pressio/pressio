
#ifndef SOLVERS_IMPL_LINE_SEARCH_MIXIN_HPP
#define SOLVERS_IMPL_LINE_SEARCH_MIXIN_HPP

#include "../../solvers_ConfigDefs.hpp"
#include "../../../../CORE_OPS"
#include "solvers_jacob_res_product_mixin.hpp"

namespace rompp{ namespace solvers{ namespace iterative{ namespace impl{

template <typename ls_tag>
struct LineSearchHelper;

template <>
struct LineSearchHelper<gn::noLineSearch>{
  template <typename scalar_t, typename ... Args>
  void operator()(scalar_t & alpha, Args&& ... args) const{
    alpha = static_cast<scalar_t>(1);
  }
};


template <>
struct LineSearchHelper<gn::ArmijoLineSearch>{

  template <typename scalar_t,   typename state_t,
	    typename residual_t, typename jacobian_t,
	    typename system_t>
  void operator()(scalar_t & alpha,
		  const state_t & y,
		  state_t & ytrial,
		  const state_t & dy,
		  residual_t & resid,
		  jacobian_t & jacob,
		  const system_t & sys) const
  {
    scalar_t c1 = 1e-4;
    alpha = static_cast<scalar_t>(1);
#ifdef DEBUG_PRINT
    ::rompp::core::io::print_stdout("line search: Armijo rule,",
				    "c1=", c1, "\n");
#endif

    ytrial.setZero();

    // eval obj function for current solution: f(y)
    auto fy  = ::rompp::core::ops::norm2(resid);

    // compute J^T * Residual
    state_t jTr(y);  jTr.setZero();
    using jtr_prod_helper_t = JacobianTranspResProdHelper<jacobian_t>;
    // evaluate
    jtr_prod_helper_t()(jacob, resid, jTr);

    // compute dy^T J^T R (this is always a dot product)
    auto c2 = ::rompp::core::ops::dot(dy, jTr);
    auto rhs = c1 * alpha * c2;

#ifdef DEBUG_PRINT
    ::rompp::core::io::print_stdout(" f(y) =", fy, "\n");
    ::rompp::core::io::print_stdout(" dy^T J^T R =", c2, "\n");
    ::rompp::core::io::print_stdout(" c1*alfa*dy^T*J^T*R =", rhs, "\n");
#endif

    bool done = false;
    while (not done)
    {
#ifdef DEBUG_PRINT
      ::rompp::core::io::print_stdout(" backtracking: alpha =",
				      alpha, "\n");
#endif

      // update
      ytrial = y + alpha * dy;

      // eval function for updated step solition: f(y + alpha*dy)
      sys.residual(ytrial, resid);
      auto fytrial  = ::rompp::core::ops::norm2(resid);
      auto lhs = fytrial-fy;

#ifdef DEBUG_PRINT
      ::rompp::core::io::print_stdout(" f(y+alpha*dy) =", fytrial, "\n");
      ::rompp::core::io::print_stdout(" f(y+alpha*dy)-f(y) =", lhs,
				      "; rhs =", rhs, "\n");
#endif

      // eval Armijo
      if (lhs <= rhs){
#ifdef DEBUG_PRINT
	::rompp::core::io::print_stdout(" lsearch done","\n");
#endif
	done = true;
      }

      // exit also when abs(fytrail-fy) < eps, leave eps = 1e-14 for now
      // change later with some machine epsilon
      if (std::abs(lhs) <= 1e-14){
#ifdef DEBUG_PRINT
	::rompp::core::io::print_stdout(" detected negligible",
					"change in obj f:",
					"abs(fytrail-fy) < 1e-14,",
					"exiting linsearch","\n");
#endif
	done = true;
      }

      /* convectional way to backtrack
       * this is equivalent to using beta^m instead of alpha
       * where m=0,1,2,...
       * and stopping when we find the smallest integer
       * to satisfy criterion
       */
      if (!done) alpha *= 0.5;

    }//while

#ifdef DEBUG_PRINT
    ::rompp::core::io::print_stdout("after line search:",
				    "alpha =", alpha, "\n");
#endif
  }//()

};


}}}} //end namespace rompp::solvers::iterative::impl
#endif






// private:
//   /* evaluating J^T*resid can be done smartly
//    * depending on the type of J.
//    * J is a matrix wrapper: then do regular mat-vec
//    * J is a multivec wrapper: use dot product
//    */
//   template <typename resid_t,
// 	    typename J_type,
// 	    typename result_t,
// 	    core::meta::enable_if_t<
// 	      core::meta::is_core_matrix_wrapper<J_type>::value
// 	      > * = nullptr>
//   void jacobTprodResidual(const resid_t & R,
// 			  const J_type & J,
// 			  result_t & result) const{
//     constexpr bool transposeJ = true;
//     ::rompp::core::ops::product<J_type, resid_t, result_t,
// 				transposeJ>(J, R, result);
//   }

//   template <typename resid_t,
// 	    typename J_type,
// 	    typename result_t,
// 	    core::meta::enable_if_t<
// 	      core::meta::is_core_multi_vector_wrapper<J_type>::value
// 	      > * = nullptr>
//   void jacobTprodResidual(const resid_t & R,
// 			  const J_type & J,
// 			  result_t & result) const{
//     ::rompp::core::ops::dot(J, R, result);
//   }
