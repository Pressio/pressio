/*
//@HEADER
// ************************************************************************
//
// solvers_create_public_api.hpp
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef SOLVERS_NONLINEAR_SOLVERS_CREATE_GAUSS_NEWTON_HPP_
#define SOLVERS_NONLINEAR_SOLVERS_CREATE_GAUSS_NEWTON_HPP_

#include "solver_impl.hpp"

namespace pressio{

namespace nonlinearsolvers{

template<class T, class = void>
struct NormalEqsDefaultTypes{
  using hessian_type  = void;
  using gradient_type = void;
};

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<class T>
struct NormalEqsDefaultTypes<
  T, mpl::enable_if_t<::pressio::is_vector_eigen<T>::value> >
{
  using hessian_type = Eigen::Matrix<typename Traits<T>::scalar_type, -1, -1>;
  using gradient_type = T;

  static hessian_type createHessian(const T & v){
    const auto ext = ::pressio::ops::extent(v, 0);
    return hessian_type(ext, ext);
  }
};
#endif

} //end nonlinearsolvers


/*
  let r(x) be the residuals with r \in R^n and x \in R^k, with k < n,
  then gauss-newton minimizes the sum of squares: S(x) = (1/2) * \sum r_i(x)*r_i(x)
  which corresponds to solving the normal equations: H delta = -g
  where H = J^T_r*J_r, and g = J^T_r * r, delta = x_k+1 - x_k.
  Here the user computes the residual and jacobian which we can use to compute
  the hssian and gradient.
*/
template<class SystemType, class LinearSolverType>
#ifdef PRESSIO_ENABLE_CXX20
requires
   (nonlinearsolvers::OverdeterminedRealValuedSystemWithResidualAndJacobian<SystemType>
 || nonlinearsolvers::OverdeterminedRealValuedSystemWithFusedResidualAndJacobian<SystemType>)
  && requires(typename SystemType::state_type & x,
	      typename SystemType::state_type & b,
	      typename SystemType::state_type & c,
	      typename SystemType::jacobian_type & J,
	      typename SystemType::residual_type & r,
	      typename nonlinearsolvers::NormalEqsDefaultTypes<typename SystemType::state_type>::hessian_type & H,
	      typename nonlinearsolvers::NormalEqsDefaultTypes<typename SystemType::state_type>::gradient_type & g,
	      nonlinearsolvers::scalar_of_t<SystemType> alpha)
  {
    { ::pressio::ops::norm2(x) } -> std::same_as< nonlinearsolvers::scalar_of_t<SystemType> >;
    { ::pressio::ops::norm2(r) } -> std::same_as< nonlinearsolvers::scalar_of_t<SystemType> >;
    { ::pressio::ops::deep_copy(b, x) };
    { ::pressio::ops::scale (x, nonlinearsolvers::scalar_of_t<SystemType>{}) };
    { ::pressio::ops::update(x, nonlinearsolvers::scalar_of_t<SystemType>{},
			     b, nonlinearsolvers::scalar_of_t<SystemType>{}) };
    { ::pressio::ops::update(x, nonlinearsolvers::scalar_of_t<SystemType>{},
			     b, nonlinearsolvers::scalar_of_t<SystemType>{},
			     c, nonlinearsolvers::scalar_of_t<SystemType>{}) };
    // to compute H and g
    { ::pressio::ops::product(transpose(), nontranspose(), 1, J, 0, H) };
    { ::pressio::ops::product(transpose(), 1, J, r, 0, g) };
  }
  && requires(typename SystemType::state_type & x,
	      typename nonlinearsolvers::NormalEqsDefaultTypes<typename SystemType::state_type>::hessian_type & H,
	      typename nonlinearsolvers::NormalEqsDefaultTypes<typename SystemType::state_type>::gradient_type & g,
	      LinearSolverType & linSolver)
  {
    linSolver.solve(H, g, x);
  }
#endif

auto create_gauss_newton_solver(const SystemType & system,
				LinearSolverType && linSolver)
{

  using state_t    = typename SystemType::state_type;
  using r_t        = typename SystemType::residual_type;
  using j_t        = typename SystemType::jacobian_type;
  using hg_default = nonlinearsolvers::NormalEqsDefaultTypes<typename SystemType::state_type>;
  using hessian_t  = typename hg_default::hessian_type;
  using gradient_t = typename hg_default::gradient_type;
  using scalar_t   = scalar_trait_t<state_t>;

  using tags = std::tuple<nonlinearsolvers::CorrectionTag,
			  nonlinearsolvers::InitialGuessTag,
			  nonlinearsolvers::ResidualTag,
			  nonlinearsolvers::JacobianTag,
			  nonlinearsolvers::GradientTag,
			  nonlinearsolvers::HessianTag,
			  nonlinearsolvers::InnerSolverTag>;
  using types = std::tuple<state_t, state_t, r_t, j_t, gradient_t, hessian_t,
			   utils::InstanceOrReferenceWrapper<LinearSolverType> >;
  using registry_t = nonlinearsolvers::impl::TagBasedStaticRegistry<tags, types>;

  registry_t reg(system.createState(), system.createState(),
		 system.createResidual(), system.createJacobian(),
		 gradient_t(system.createState()),
		 hessian_t( hg_default::createHessian(system.createState()) ),
		 std::forward<LinearSolverType>(linSolver));

  using nonlinearsolvers::Diagnostic;
  const std::vector<Diagnostic> defaultDiagnostics =
    {Diagnostic::objectiveAbsolute,
     Diagnostic::objectiveRelative,
     Diagnostic::correctionAbsolutel2Norm,
     Diagnostic::correctionRelativel2Norm,
     Diagnostic::gradientAbsolutel2Norm,
     Diagnostic::gradientRelativel2Norm};

  using tag = nonlinearsolvers::impl::GaussNewtonNormalEqTag;
  return nonlinearsolvers::impl::NonLinLeastSquares<tag, state_t, registry_t, scalar_t>
    (tag{}, std::move(reg), defaultDiagnostics);
}


/*
  let r(x) be the residuals with r \in R^n and x \in R^k, with k < n,
  then weighted gauss-newton minimizes the weighted sum of squares:
      S(x) = (1/2) * \sum r_i(x) * w_i * r_i(x)
  which corresponds to solving the normal equations: H delta = -g
  where:
     delta = x_k+1 - x_k
  and
     H = J^T_r* W * J_r
     g = J^T_r * W * r
*/
template<class SystemType, class LinearSolverType, class WeightingOpType>
#ifdef PRESSIO_ENABLE_CXX20
requires
     (nonlinearsolvers::OverdeterminedRealValuedSystemWithResidualAndJacobian<SystemType>
   || nonlinearsolvers::OverdeterminedRealValuedSystemWithFusedResidualAndJacobian<SystemType>)
  && requires(typename SystemType::state_type & x,
	      typename SystemType::state_type & b,
	      typename SystemType::state_type & c,
	      const typename SystemType::jacobian_type & J,
	      const typename SystemType::residual_type & r,
	      typename SystemType::jacobian_type & WJ,
	      typename SystemType::residual_type & Wr,
	      typename nonlinearsolvers::NormalEqsDefaultTypes<typename SystemType::state_type>::hessian_type & H,
	      typename nonlinearsolvers::NormalEqsDefaultTypes<typename SystemType::state_type>::gradient_type & g,
	      nonlinearsolvers::scalar_of_t<SystemType> alpha,
	      const WeightingOpType & W)
  {
    { ::pressio::ops::norm2(x) } -> std::same_as< nonlinearsolvers::scalar_of_t<SystemType> >;
    { ::pressio::ops::norm2(r) } -> std::same_as< nonlinearsolvers::scalar_of_t<SystemType> >;
    { ::pressio::ops::deep_copy(b, x) };
    { ::pressio::ops::scale (x, nonlinearsolvers::scalar_of_t<SystemType>{}) };
    { ::pressio::ops::update(x, nonlinearsolvers::scalar_of_t<SystemType>{},
			     b, nonlinearsolvers::scalar_of_t<SystemType>{}) };
    { ::pressio::ops::update(x, nonlinearsolvers::scalar_of_t<SystemType>{},
			     b, nonlinearsolvers::scalar_of_t<SystemType>{},
			     c, nonlinearsolvers::scalar_of_t<SystemType>{}) };
    { ::pressio::ops::dot(r, Wr) } -> std::same_as< nonlinearsolvers::scalar_of_t<SystemType> >;

    // how we use the weighting
    { W(r, Wr) };
    { W(J, WJ) };

    // to compute H and g
    { ::pressio::ops::product(transpose(), nontranspose(), 1, J, WJ, 0, H) };
    { ::pressio::ops::product(transpose(), 1, J, Wr, 0, g) };
  }
  && requires(typename SystemType::state_type & x,
	      typename nonlinearsolvers::NormalEqsDefaultTypes<typename SystemType::state_type>::hessian_type & H,
	      typename nonlinearsolvers::NormalEqsDefaultTypes<typename SystemType::state_type>::gradient_type & g,
	      LinearSolverType & linSolver)
  {
    linSolver.solve(H, g, x);
  }
#endif

auto create_gauss_newton_solver(const SystemType & system,
				LinearSolverType && linSolver,
				WeightingOpType && W)
{

  using state_t    = typename SystemType::state_type;
  using r_t        = typename SystemType::residual_type;
  using j_t        = typename SystemType::jacobian_type;
  using hg_default = nonlinearsolvers::NormalEqsDefaultTypes<typename SystemType::state_type>;
  using hessian_t  = typename hg_default::hessian_type;
  using gradient_t = typename hg_default::gradient_type;
  using scalar_t   = scalar_trait_t<state_t>;

  using tags = std::tuple<nonlinearsolvers::CorrectionTag,
			  nonlinearsolvers::InitialGuessTag,
			  nonlinearsolvers::ResidualTag,
			  nonlinearsolvers::JacobianTag,
			  nonlinearsolvers::WeightedResidualTag,
			  nonlinearsolvers::WeightedJacobianTag,
			  nonlinearsolvers::GradientTag,
			  nonlinearsolvers::HessianTag,
			  nonlinearsolvers::InnerSolverTag,
			  nonlinearsolvers::WeightingOperatorTag>;
  using types = std::tuple<state_t, state_t,
			   r_t, j_t,
			   r_t, j_t, // weighted r and J have same type of r,J
			   gradient_t, hessian_t,
			   utils::InstanceOrReferenceWrapper<LinearSolverType>,
			   utils::InstanceOrReferenceWrapper<WeightingOpType>>;
  using registry_t = nonlinearsolvers::impl::TagBasedStaticRegistry<tags, types>;

  registry_t reg(system.createState(), system.createState(),
		 system.createResidual(), system.createJacobian(),
		 system.createResidual(), system.createJacobian(),
		 gradient_t(system.createState()),
		 hessian_t( hg_default::createHessian(system.createState()) ),
		 std::forward<LinearSolverType>(linSolver),
		 std::forward<WeightingOpType>(W));

  using nonlinearsolvers::Diagnostic;
  const std::vector<Diagnostic> defaultDiagnostics =
    {Diagnostic::objectiveAbsolute,
     Diagnostic::objectiveRelative,
     Diagnostic::correctionAbsolutel2Norm,
     Diagnostic::correctionRelativel2Norm,
     Diagnostic::gradientAbsolutel2Norm,
     Diagnostic::gradientRelativel2Norm};

  using tag = nonlinearsolvers::impl::WeightedGaussNewtonNormalEqTag;
  return nonlinearsolvers::impl::NonLinLeastSquares<tag, state_t, registry_t, scalar_t>
    (tag{}, std::move(reg), defaultDiagnostics);
}


} // end namespace pressio
#endif  // SOLVERS_NONLINEAR_SOLVERS_CREATE_PUBLIC_API_HPP_








// template<class SystemType, class LinearSolverType>
// #ifdef PRESSIO_ENABLE_CXX20
// requires
//    (SystemWithHessianAndGradient<SystemType>
//    || SystemWithFusedHessianAndGradient<SystemType>)
//   && LinearSolver<
//       mpl::remove_cvref_t<LinearSolverType>,
//       typename SystemType::state_type,
//       typename SystemType::hessian_type,
//       typename SystemType::gradient_type
//      >
// #endif
// auto create_gauss_newton(const SystemType & system,
// 			 LinearSolverType && linSolver)
// {
//   using linear_solver_type = mpl::remove_cvref_t<LinearSolverType>;

//   using system_type = SystemType;
//   using state_t    = typename system_type::state_type;
//   using hessian_t   = typename system_type::hessian_type;
//   using gradient_t  = typename system_type::gradient_type;
//   using norm_value_type = decltype(::pressio::ops::norm2(std::declval<gradient_t>()));

//   using tags = std::tuple<CorrectionTag, StateTag,
// 			  InitialGuessTag, LineSearchTrialStateTag,
// 			  GradientTag, HessianTag,
// 			  InnerSolverTag>;
//   using types = std::tuple<state_t, state_t,
// 			   state_t, state_t,
// 			   gradient_t, hessian_t,
// 			   utils::InstanceOrReferenceWrapper<LinearSolverType>>;
//   using registry_t = RegistryImpl<tags, types>;
//   registry_t reg(system.createState(), system.createState(),
// 		 system.createState(), system.createState(),
// 		 system.createGradient(), system.createHessian(),
// 		 std::forward<LinearSolverType>(linSolver));

//   const std::vector<Diagnostic> defaultDiagnostics =
//     {Diagnostic::residualAbsolutel2Norm,
//      Diagnostic::residualRelativel2Norm,
//      Diagnostic::correctionAbsolutel2Norm,
//      Diagnostic::correctionRelativel2Norm,
//      Diagnostic::gradientAbsolutel2Norm,
//      Diagnostic::gradientRelativel2Norm};

//   return NonLinLs<
//     GaussNewtonNormalEqTag, state_t,
//     registry_t, norm_value_type>(GaussNewtonNormalEqTag{}, reg, defaultDiagnostics);
// }


// template<class SystemType, class QRSolverType>
// #ifdef PRESSIO_ENABLE_CXX20
// requires
//      (OverdeterminedSystemWithResidualAndJacobian<SystemType>
//    || OverdeterminedSystemWithFusedResidualAndJacobian<SystemType>)
//   && QRSolverForGnQr<
// 	mpl::remove_cvref_t<QRSolverType>,
// 	typename SystemType::state_type,
// 	typename SystemType::jacobian_type,
// 	typename SystemType::residual_type
//      >
// #endif
// auto create_gauss_newtonQR(const SystemType & system,
// 			   QRSolverType && linSolver)
// {
//   using system_type = SystemType;
//   using state_t    = typename system_type::state_type;
//   using r_t        = typename system_type::residual_type;
//   using j_t        = typename system_type::jacobian_type;
//   using gradient_t = state_t;
//   using norm_value_type = decltype(::pressio::ops::norm2(std::declval<gradient_t>()));

//   using tags = std::tuple<CorrectionTag, StateTag,
// 			  InitialGuessTag, LineSearchTrialStateTag,
// 			  ResidualTag, JacobianTag,
// 			  GradientTag,
// 			  InnerSolverTag>;
//   using types = std::tuple<state_t, state_t,
// 			   state_t, state_t,
// 			   r_t, j_t,
// 			   gradient_t,
// 			   utils::InstanceOrReferenceWrapper<QRSolverType>>;
//   using registry_t = RegistryImpl<tags, types>;
//   registry_t reg(system.createState(), system.createState(),
// 		 system.createState(), system.createState(),
// 		 system.createResidual(), system.createJacobian(),
// 		 gradient_t(system.createState()),
// 		 std::forward<QRSolverType>(linSolver));

//   const std::vector<Diagnostic> defaultDiagnostics =
//     {Diagnostic::residualAbsolutel2Norm,
//      Diagnostic::residualRelativel2Norm,
//      Diagnostic::correctionAbsolutel2Norm,
//      Diagnostic::correctionRelativel2Norm,
//      Diagnostic::gradientAbsolutel2Norm,
//      Diagnostic::gradientRelativel2Norm};

//   return NonLinLs<
//     GaussNewtonQRTag, state_t,
//     registry_t, norm_value_type>(GaussNewtonQRTag{}, reg, defaultDiagnostics);
// }


// // ----------------------------------------------------------------
// // Leven-Marq
// // ----------------------------------------------------------------

// template<class SystemType, class LinearSolverType>
// #ifdef PRESSIO_ENABLE_CXX20
// requires
//      (OverdeterminedSystemWithResidualAndJacobian<SystemType>
//    || OverdeterminedSystemWithFusedResidualAndJacobian<SystemType>)
//   && LinearSolver<
//       mpl::remove_cvref_t<LinearSolverType>,
//       typename SystemType::state_type,
//       typename NormalEqsDefaultTypes<typename SystemType::state_type>::hessian_type,
//       typename NormalEqsDefaultTypes<typename SystemType::state_type>::gradient_type
//      >
// #endif
// auto create_levenberg_marquardt(const SystemType & system,
// 				LinearSolverType && linSolver)
// {
//   using linear_solver_type = mpl::remove_cvref_t<LinearSolverType>;

//   using system_type = SystemType;
//   using state_t    = typename system_type::state_type;
//   using r_t        = typename system_type::residual_type;
//   using j_t        = typename system_type::jacobian_type;
//   using hg_default = NormalEqsDefaultTypes<typename SystemType::state_type>;
//   using hessian_t  = typename hg_default::hessian_type;
//   using gradient_t = typename hg_default::gradient_type;
//   using norm_value_type = decltype(::pressio::ops::norm2(std::declval<gradient_t>()));

//   using lm_damp_t = typename ::pressio::Traits<hessian_t>::scalar_type;

//   using tags = std::tuple<CorrectionTag, StateTag,
// 			  InitialGuessTag, LineSearchTrialStateTag,
// 			  ResidualTag, JacobianTag,
// 			  GradientTag, HessianTag,
// 			  LevenbergMarquardtUndampedHessianTag,
// 			  LevenbergMarquardtDampingTag,
// 			  InnerSolverTag>;
//   using types = std::tuple<state_t, state_t,
// 			   state_t, state_t,
// 			   r_t, j_t,
// 			   gradient_t, hessian_t,
// 			   hessian_t,
// 			   lm_damp_t,
// 			   utils::InstanceOrReferenceWrapper<LinearSolverType>>;
//   using registry_t = RegistryImpl<tags, types>;

//   const auto stateExt = ::pressio::ops::extent(system.createState(), 0);
//   registry_t reg(system.createState(), system.createState(),
// 		 system.createState(), system.createState(),
// 		 system.createResidual(), system.createJacobian(),
// 		 gradient_t(system.createState()),
// 		 hessian_t(stateExt, stateExt),
// 		 hessian_t(stateExt, stateExt),
// 		 lm_damp_t{1},
// 		 std::forward<LinearSolverType>(linSolver));

//   const std::vector<Diagnostic> defaultDiagnostics =
//     {Diagnostic::residualAbsolutel2Norm,
//      Diagnostic::residualRelativel2Norm,
//      Diagnostic::correctionAbsolutel2Norm,
//      Diagnostic::correctionRelativel2Norm,
//      Diagnostic::gradientAbsolutel2Norm,
//      Diagnostic::gradientRelativel2Norm};

//   return NonLinLs<
//     LevenbergMarquardtNormalEqTag, state_t,
//     registry_t, norm_value_type>(LevenbergMarquardtNormalEqTag{}, reg, defaultDiagnostics);
// }
