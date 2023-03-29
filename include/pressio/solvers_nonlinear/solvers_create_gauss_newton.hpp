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

/*
  Gauss-newton minimizes the sum of squares: S(x) = (1/2) * \sum r_i(x)*r_i(x)
  where r(x) is the residual vector with r \in R^n and x \in R^k, with k < n,
  which corresponds to solving the normal equations: H delta = -g
  where H = J^T_r*J_r, and g = J^T_r * r, delta = x_k+1 - x_k.
  The user provides a problem computing residual and jacobian
  which we can use to compute the hessian and gradient.
*/

template<class SystemType, class LinearSolverType>
#ifdef PRESSIO_ENABLE_CXX20
  requires nonlinearsolvers::RealValuedNonlinearSystemFusingResidualAndJacobian<SystemType>
  && nonlinearsolvers::valid_state_for_least_squares<typename SystemType::state_type>::value
  && (Traits<typename SystemType::state_type>::rank == 1)
  && (Traits<typename SystemType::residual_type>::rank == 1)
  && (Traits<typename SystemType::jacobian_type>::rank == 2)
  && requires(typename SystemType::state_type & x,
	      const typename SystemType::jacobian_type & J,
	      const typename SystemType::residual_type & r,
	      nonlinearsolvers::normal_eqs_default_hessian_t<typename SystemType::state_type>  & H,
	      nonlinearsolvers::normal_eqs_default_gradient_t<typename SystemType::state_type> & g,
	      LinearSolverType && linSolver)
  {
    { ::pressio::ops::norm2(r) } -> std::same_as< nonlinearsolvers::scalar_of_t<SystemType> >;
    { ::pressio::ops::product(transpose(), nontranspose(), 1, J, 0, H) };
    { ::pressio::ops::product(transpose(), 1, J, r, 0, g) };
    { linSolver.solve(std::as_const(H), std::as_const(g), x) };
  }
#endif
auto create_gauss_newton_solver(const SystemType & system,           /*(1)*/
				LinearSolverType && linSolver)
{
  using scalar_t = nonlinearsolvers::scalar_of_t<SystemType>;
  using state_t    = typename SystemType::state_type;
  using r_t        = typename SystemType::residual_type;
  using j_t        = typename SystemType::jacobian_type;
  using hg_default = nonlinearsolvers::normal_eqs_default_types<state_t>;
  using hessian_t  = typename hg_default::hessian_type;
  using gradient_t = typename hg_default::gradient_type;

  using tags = std::tuple<
    nonlinearsolvers::CorrectionTag,
    nonlinearsolvers::InitialGuessTag,
    nonlinearsolvers::ResidualTag,
    nonlinearsolvers::JacobianTag,
    nonlinearsolvers::GradientTag,
    nonlinearsolvers::HessianTag,
    nonlinearsolvers::InnerSolverTag,
    nonlinearsolvers::impl::SystemTag
    >;
  using types = std::tuple<
    state_t, state_t, r_t, j_t, gradient_t, hessian_t,
    utils::InstanceOrReferenceWrapper<LinearSolverType>,
    SystemType const *
    >;

  using registry_t = nonlinearsolvers::impl::TagBasedStaticRegistry<tags, types>;
  registry_t reg(system.createState(), system.createState(),
		 system.createResidual(), system.createJacobian(),
		 gradient_t(system.createState()),
		 hessian_t( hg_default::createHessian(system.createState()) ),
		 std::forward<LinearSolverType>(linSolver),
		 &system);

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
  Weighted gauss-newton minimizes the weighted sum of squares:
      S(x) = (1/2) * \sum r_i(x) * w_i * r_i(x)
  where r(x) is the residual vector with r \in R^n and x \in R^k, with k < n,
  which corresponds to solving the normal equations: H delta = -g
  where: delta = x_k+1 - x_k, H = J^T_r* W * J_r, and g = J^T_r * W * r
*/

template<class SystemType, class LinearSolverType, class WeightingOpType>
#ifdef PRESSIO_ENABLE_CXX20
  requires nonlinearsolvers::RealValuedNonlinearSystemFusingResidualAndJacobian<SystemType>
  && nonlinearsolvers::valid_state_for_least_squares<typename SystemType::state_type>::value
  && (Traits<typename SystemType::state_type>::rank == 1)
  && (Traits<typename SystemType::residual_type>::rank == 1)
  && (Traits<typename SystemType::jacobian_type>::rank == 2)
  && requires(typename SystemType::state_type & x,
	      const typename SystemType::jacobian_type & J,
	      const typename SystemType::residual_type & r,
	      typename SystemType::jacobian_type & WJ,
	      typename SystemType::residual_type & Wr,
	      nonlinearsolvers::normal_eqs_default_hessian_t<typename SystemType::state_type>  & H,
	      nonlinearsolvers::normal_eqs_default_gradient_t<typename SystemType::state_type> & g,
	      WeightingOpType && W,
	      LinearSolverType && linSolver)
  {
    { ::pressio::ops::norm2(r) }
	-> std::same_as< nonlinearsolvers::scalar_of_t<SystemType> >;
    { ::pressio::ops::dot(std::as_const(r), std::as_const(Wr)) }
	-> std::same_as< nonlinearsolvers::scalar_of_t<SystemType> >;

    { W(r, Wr) };
    { W(J, WJ) };
    { ::pressio::ops::product(transpose(), nontranspose(), 1, J, std::as_const(WJ), 0, H) };
    { ::pressio::ops::product(transpose(), 1, J, std::as_const(Wr), 0, g) };

    { linSolver.solve(std::as_const(H), std::as_const(g), x) };
  }
#endif
auto create_gauss_newton_solver(const SystemType & system,           /*(2)*/
				LinearSolverType && linSolver,
				WeightingOpType && weigher)
{
  using scalar_t   = nonlinearsolvers::scalar_of_t<SystemType>;
  using state_t    = typename SystemType::state_type;
  using r_t        = typename SystemType::residual_type;
  using j_t        = typename SystemType::jacobian_type;
  using hg_default = nonlinearsolvers::normal_eqs_default_types<state_t>;
  using hessian_t  = typename hg_default::hessian_type;
  using gradient_t = typename hg_default::gradient_type;

  using tags = std::tuple<
    nonlinearsolvers::CorrectionTag,
    nonlinearsolvers::InitialGuessTag,
    nonlinearsolvers::ResidualTag,
    nonlinearsolvers::JacobianTag,
    nonlinearsolvers::WeightedResidualTag,
    nonlinearsolvers::WeightedJacobianTag,
    nonlinearsolvers::GradientTag,
    nonlinearsolvers::HessianTag,
    nonlinearsolvers::InnerSolverTag,
    nonlinearsolvers::WeightingOperatorTag,
    nonlinearsolvers::impl::SystemTag
    >;
  using types = std::tuple<
    state_t, state_t,
    r_t, j_t, r_t, j_t, // Wr and WJ have same type of r,J
    gradient_t, hessian_t,
    utils::InstanceOrReferenceWrapper<LinearSolverType>,
    utils::InstanceOrReferenceWrapper<WeightingOpType>,
    SystemType const *
    >;
  using registry_t = nonlinearsolvers::impl::TagBasedStaticRegistry<tags, types>;

  registry_t reg(system.createState(), system.createState(),
		 system.createResidual(), system.createJacobian(),
		 system.createResidual(), system.createJacobian(),
		 gradient_t(system.createState()),
		 hessian_t( hg_default::createHessian(system.createState()) ),
		 std::forward<LinearSolverType>(linSolver),
		 std::forward<WeightingOpType>(weigher),
		 &system);

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

namespace experimental{

/*
  see: https://en.wikipedia.org/wiki/Non-linear_least_squares#QR_decomposition
  but careful that here we use J to refer to jacobian wrt r,
  which is the negative of J_f.

  NOTE: this is inside the experimental namespace because I know this works as is
  also for Trilinos, but I am not happy with the qr solver part, I think that can be
  improved but we can do that later since it needs to be synced up with the refactoring
  of the pressio/qr component
*/

template<class SystemType, class QRSolverType>
#ifdef PRESSIO_ENABLE_CXX20
  requires nonlinearsolvers::RealValuedNonlinearSystemFusingResidualAndJacobian<SystemType>
  && nonlinearsolvers::valid_state_for_least_squares< typename SystemType::state_type >::value
  && (Traits<typename SystemType::state_type>::rank    == 1)
  && (Traits<typename SystemType::residual_type>::rank == 1)
  && (Traits<typename SystemType::jacobian_type>::rank == 2)
  && requires(typename SystemType::state_type & x,
	      const typename SystemType::jacobian_type & J,
	      const typename SystemType::residual_type & r,
	      typename SystemType::state_type & QTr,
	      QRSolverType && qrSolver)
  {
    { ::pressio::ops::norm2(r) } -> std::same_as< nonlinearsolvers::scalar_of_t<SystemType> >;
    { qrSolver.computeThin(J)               } -> std::same_as<void>;
    { qrSolver.applyQTranspose(r, QTr)      } -> std::same_as<void>;
    { qrSolver.solve(std::as_const(QTr), x) } -> std::same_as<void>;
    //{ A.applyRTranspose(d, b) } -> std::same_as<void>;
  }
#endif
auto create_gauss_newton_qr_solver(const SystemType & system,
				   QRSolverType && linSolver)
{
  using scalar_t = nonlinearsolvers::scalar_of_t<SystemType>;
  using state_t    = typename SystemType::state_type;
  using r_t        = typename SystemType::residual_type;
  using j_t        = typename SystemType::jacobian_type;
  using QTr_t      = state_t; // type of Q^T*r
  using gradient_t = state_t; // type of J^T r

  using tags = std::tuple<
    nonlinearsolvers::CorrectionTag,
    nonlinearsolvers::InitialGuessTag,
    nonlinearsolvers::ResidualTag,
    nonlinearsolvers::JacobianTag,
    nonlinearsolvers::GradientTag,
    nonlinearsolvers::impl::QTransposeResidualTag,
    nonlinearsolvers::InnerSolverTag,
    nonlinearsolvers::impl::SystemTag
    >;
  using types = std::tuple<
    state_t, state_t, r_t, j_t, gradient_t, QTr_t,
    utils::InstanceOrReferenceWrapper<QRSolverType>,
    SystemType const *
    >;

  using registry_t = nonlinearsolvers::impl::TagBasedStaticRegistry<tags, types>;
  registry_t reg(system.createState(), system.createState(),
		 system.createResidual(), system.createJacobian(),
		 system.createState(), // gradient is same extent as state
		 system.createState(), //Q^T*r has same extent as state
		 std::forward<QRSolverType>(linSolver),
		 &system);

  using nonlinearsolvers::Diagnostic;
  const std::vector<Diagnostic> defaultDiagnostics =
    {Diagnostic::objectiveAbsolute,
     Diagnostic::objectiveRelative,
     Diagnostic::correctionAbsolutel2Norm,
     Diagnostic::correctionRelativel2Norm,
     Diagnostic::gradientAbsolutel2Norm,
     Diagnostic::gradientRelativel2Norm};

  using tag = nonlinearsolvers::impl::GaussNewtonQrTag;
  return nonlinearsolvers::impl::NonLinLeastSquares<tag, state_t, registry_t, scalar_t>
    (tag{}, std::move(reg), defaultDiagnostics);
}
}//end experimental


} // end namespace pressio
#endif  // SOLVERS_NONLINEAR_SOLVERS_CREATE_PUBLIC_API_HPP_





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
