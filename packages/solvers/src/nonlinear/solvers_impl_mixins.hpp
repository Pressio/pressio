
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
  bool operator()(const state_t & x, const state_t & dx,
		  scalar_t norm_dx, step_t step,
		  step_t maxIters, scalar_t tol) const{
    return step==maxIters;
  }
};

template <typename norm_t>
struct isConvergedHelper<
  converged_when::absoluteNormCorrectionBelowTol<norm_t>>{

  static constexpr char const * description_ = "norm(dx) < tol";

  template <typename state_t, typename step_t, typename scalar_t>
  bool operator()(const state_t & x, const state_t & dx,
		  scalar_t norm_dx, step_t step,
		  step_t maxIters, scalar_t tol) const{
    return (norm_dx<tol);
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
		  const state_t & x,
		  const state_t & dx,
		  residual_t & resid,
		  jacobian_t & jacob,
		  const system_t & sys) const
  {
    scalar_t c1 = 1e-4;
    alpha = static_cast<scalar_t>(1);
#ifdef DEBUG_PRINT
    ::rompp::core::io::print_stdout("line search: Armijo rule", "\n");
#endif

    // create a new state vector for doing line search
    state_t xtrial(x);

    // compute J^T * Residual
    state_t jTr(x);
    jTr.setZero();
    this->jacobTprodResidual(resid, jacob, jTr);
    //::rompp::core::ops::dot(jacob, resid, jTr);
    //*jTr.data() = (*jacob.data()).transpose() * (*resid.data());

    // compute dx^T J^T R
    auto c2 = ::rompp::core::ops::dot(dx, jTr);
    ::rompp::core::io::print_stdout(" dx^T J^T R =", c2, "\n");
    auto rhs = c1 * alpha * c2;
    ::rompp::core::io::print_stdout(" c1*alfa*dx^T*J^T*R =", rhs, "\n");

    // eval obj function for current solution: f(x)
    auto fx  = ::rompp::core::ops::norm2(resid);
    ::rompp::core::io::print_stdout(" f(x) =", fx, "\n");

    bool done = false;
    while (not done){
      ::rompp::core::io::print_stdout(" backtracking: alpha =",
				      alpha, "\n");
      // update
      xtrial = x + alpha * dx;

      // eval function for updated step solition: f(x + alpha*dx)
      sys.residual(xtrial, resid);
      auto fxtrial  = ::rompp::core::ops::norm2(resid);
      ::rompp::core::io::print_stdout(" f(x+alpha*dx) =", fxtrial, "\n");

      auto lhs = fxtrial-fx;

      ::rompp::core::io::print_stdout(" f(x+alpha*dx)-f(x) =", lhs,
				      " rhs =", rhs, "\n");
      // eval Armijo
      if (lhs <= rhs){
	::rompp::core::io::print_stdout(" lsearch done","\n");
	done = true;
      }

      // exit also when abs(fxtrail-fx) < eps, leave eps = 1e-14 for now
      // change later with some machine epsilon
      if (std::abs(lhs) <= 1e-14)
	done = true;

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
