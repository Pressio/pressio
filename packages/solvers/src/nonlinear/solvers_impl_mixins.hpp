
#ifndef SOLVERS_IMPL_MIXINS_HPP
#define SOLVERS_IMPL_MIXINS_HPP

#include "../solvers_ConfigDefs.hpp"
#include "../../../CORE_OPS"

namespace rompp{ namespace solvers{ namespace iterative{ namespace impl{

template <typename convergence_tag>
struct norm_type_to_use{
  using norm_t = L2Norm;
};

template <typename norm_type>
struct norm_type_to_use<
  converged_when::absoluteNormCorrectionBelowTol<norm_type>
  >{
  using norm_t = norm_type;
};
//---------------------------------------------------------

template <typename norm_t>
struct computeNormHelper;

template <>
struct computeNormHelper<::rompp::solvers::L2Norm>{
  template <typename vec_t, typename scalar_t>
  void operator()(const vec_t & vecIn, scalar_t & result) const{
    result = ::rompp::core::ops::norm2(vecIn);
  }
};

template <>
struct computeNormHelper<::rompp::solvers::L1Norm>{
  template <typename vec_t, typename scalar_t>
  void operator()(const vec_t & vecIn, scalar_t & result) const{
    result = ::rompp::core::ops::norm1(vecIn);
  }
};
//---------------------------------------------------------


template <typename conv_tag>
struct isConvergedHelper;

template <>
struct isConvergedHelper<converged_when::completingNumMaxIters>{
  static constexpr char const * description_ = "complete max iters";

  template <typename state_t, typename step_t, typename scalar_t>
  bool operator()(const state_t & y, const state_t & dy,
		  scalar_t norm_dy, step_t step,
		  step_t maxIters, scalar_t tol) const{
    return step==maxIters;
  }
};

template <typename norm_t>
struct isConvergedHelper<
  converged_when::absoluteNormCorrectionBelowTol<norm_t>>{

  static constexpr char const * description_ = "norm(dx) < tol";

  template <typename state_t, typename step_t, typename scalar_t>
  bool operator()(const state_t & y, const state_t & dy,
		  scalar_t norm_dy, step_t step,
		  step_t maxIters, scalar_t tol) const{
    return (norm_dy<tol);
  }
};
//---------------------------------------------------------



template <typename ls_tag>
struct lineSearchHelper;


template <>
struct lineSearchHelper<gn::noLineSearch>{
  template <typename scalar_t, typename ... Args>
  void operator()(scalar_t & alpha, Args&& ... args) const{
    alpha = static_cast<scalar_t>(1);
  }
};


template <>
struct lineSearchHelper<gn::ArmijoLineSearch>{

private:
  /* evaluating J^T*resid can be done smartly
   * depending on the type of J.
   * J is a matrix wrapper: then do regular mat-vec
   * J is a multivec wrapper: use dot product
   */
  template <typename resid_t,
	    typename J_type,
	    typename result_t,
	    core::meta::enable_if_t<
	      core::meta::is_core_matrix_wrapper<J_type>::value
	      > * = nullptr>
  void jacobTprodResidual(const resid_t & R,
			  const J_type & J,
			  result_t & result) const{
    constexpr bool transposeJ = true;
    ::rompp::core::ops::product<J_type, resid_t, result_t,
				transposeJ>(J, R, result);
  }

  template <typename resid_t,
	    typename J_type,
	    typename result_t,
	    core::meta::enable_if_t<
	      core::meta::is_core_multi_vector_wrapper<J_type>::value
	      > * = nullptr>
  void jacobTprodResidual(const resid_t & R,
			  const J_type & J,
			  result_t & result) const{
    ::rompp::core::ops::dot(J, R, result);
  }

public:
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
    state_t jTr(y);
    jTr.setZero();
    this->jacobTprodResidual(resid, jacob, jTr);

    // compute dy^T J^T R
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
				    "alpha =", alpha,
				    "\n");
#endif
  }//()

};


}}}} //end namespace rompp::solvers::iterative::impl
#endif
