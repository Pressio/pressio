
#ifndef ROM_GALERKIN_TYPE_GENERATORS_HPP_
#define ROM_GALERKIN_TYPE_GENERATORS_HPP_

#include "../rom_ConfigDefs.hpp"
#include "../rom_forward_declarations.hpp"
#include "../rom_data_fom_rhs.hpp"
#include "../rom_data_fom_states.hpp"
#include "../policies/rom_evaluate_fom_rhs_policy.hpp"
#include "../policies/rom_apply_fom_jacobian_policy.hpp"
#include "../../../ode/src/ode_forward_declarations.hpp"

namespace rompp{ namespace rom{

template <typename fom_type,
	  typename decoder_type,
	  typename galerkin_state_type,
	  typename galerkin_residual_type>
struct GalerkinCommonTypes{
  // these are native types of the full-order model (fom)
  using fom_t			= fom_type;
  using scalar_t		= typename fom_t::scalar_type;
  using fom_state_t		= typename fom_t::state_type;
  using fom_rhs_t		= typename fom_t::residual_type;

  // declare fom wrapper types
  using fom_state_w_t		= ::rompp::core::Vector<fom_state_t>;
  using fom_rhs_w_t		= ::rompp::core::Vector<fom_rhs_t>;

  // decoder types (passed in)
  using decoder_t		= decoder_type;
  using decoder_jac_t		= typename decoder_t::jacobian_t;

  // rom state type (passed in)
  using galerkin_state_t	= galerkin_state_type;

  // the GALERKIN residual type (passed in)
  using galerkin_residual_t	= galerkin_residual_type;

  // class type holding fom states data
  using fom_states_data = ::rompp::rom::FomStatesData<
	fom_state_w_t, 0, decoder_t>;

  // class type holding fom rhs data
  using fom_rhs_data = ::rompp::rom::FomRhsData<fom_rhs_w_t>;
};


template <typename fom_type,
	  ode::ExplicitEnum odeName,
	  typename decoder_type,
	  typename galerkin_state_type,
	  typename galerkin_residual_type = galerkin_state_type>
struct DefaultGalerkinTypeGenerator
  : GalerkinCommonTypes<fom_type, decoder_type,
			galerkin_state_type, galerkin_residual_type>{

  using base_t = GalerkinCommonTypes<fom_type, decoder_type,
				     galerkin_state_type,
				     galerkin_residual_type>;

  using typename base_t::fom_t;
  using typename base_t::scalar_t;
  using typename base_t::fom_state_t;
  using typename base_t::fom_state_w_t;
  using typename base_t::fom_rhs_w_t;
  using typename base_t::decoder_t;
  using typename base_t::decoder_jac_t;
  using typename base_t::galerkin_state_t;
  using typename base_t::galerkin_residual_t;
  using typename base_t::fom_states_data;
  using typename base_t::fom_rhs_data;

  // policy for evaluating the ode residual
  using galerkin_residual_policy_t =
    ::rompp::rom::DefaultGalerkinExplicitResidualPolicy<
    fom_states_data, fom_rhs_data, decoder_t>;

  // declare type of stepper object
  using galerkin_stepper_t = ::rompp::ode::ExplicitStepper<
    odeName, galerkin_state_type, fom_type,
    galerkin_residual_t, galerkin_residual_policy_t>;

};//end class





// template <typename fom_type,
// 	  ode::ImplicitEnum odeName,
// 	  typename decoder_type,
// 	  typename lspg_state_type>
// struct PreconditionedLSPGTypeGenerator
//   : LSPGCommonTypes<fom_type, odeName, decoder_type, lspg_state_type>{

//   using base_t = LSPGCommonTypes<fom_type, odeName, decoder_type, lspg_state_type>;

//   using typename base_t::fom_t;
//   using typename base_t::scalar_t;
//   using typename base_t::fom_state_t;
//   using typename base_t::fom_state_w_t;
//   using typename base_t::fom_rhs_w_t;
//   using typename base_t::decoder_t;
//   using typename base_t::decoder_jac_t;
//   using typename base_t::lspg_state_t;
//   using typename base_t::lspg_residual_t;
//   using typename base_t::fom_states_data;
//   using typename base_t::fom_rhs_data;

//   /* lspg_matrix_t is type of J*decoder_jac_t (in the most basic case) where
//    * * J is the jacobian of the fom rhs
//    * * decoder_jac_t is the type of the decoder jacobian
//    * In more complex cases, we might have (something)*J*decoder_jac_t,
//    * where (something) is product of few matrices.
//    * For now, set lspg_matrix_t to be of same type as decoder_jac_t
//    * if phi is MV<>, then lspg_matrix_t = core::MV<>
//    * if phi is Matrix<>, then we have core::Matrix<>
//    * not a bad assumption since all matrices are left-applied to decoder_jac_t
//    */
//   using lspg_matrix_t		= decoder_jac_t;

//   // policy for evaluating the rhs of the fom object
//   using fom_eval_rhs_policy_t	= rom::policy::EvaluateFomRhsDefault;

//   // policy for left multiplying the fom jacobian with decoder_jac_t
//   // possibly involving other stuff like explained above
//   using fom_apply_jac_policy_t	= rom::policy::ApplyFomJacobianDefault;

//   // policy defining how to compute the LSPG time-discrete residual
//   using lspg_residual_policy_t =
//     rom::decorator::Preconditioned<
//     rom::LSPGResidualPolicy<
//       fom_states_data, fom_rhs_data, fom_eval_rhs_policy_t
//       >
//     >;

//   // policy defining how to compute the LSPG time-discrete jacobian
//   using lspg_jacobian_policy_t	=
//     rom::decorator::Preconditioned<
//     rom::LSPGJacobianPolicy<
//       fom_states_data, lspg_matrix_t, fom_apply_jac_policy_t
//       >
//     >;

//   using aux_stepper_t = typename auxStepperHelper<odeName, lspg_state_type,
// 						  lspg_residual_t, lspg_matrix_t,
// 						  fom_type, lspg_residual_policy_t,
// 						  lspg_jacobian_policy_t>::type;

//   // declare type of stepper object
//   using rom_stepper_t		= ode::ImplicitStepper<
//     odeName, lspg_state_type, lspg_residual_t, lspg_matrix_t,
//     fom_type, aux_stepper_t, lspg_residual_policy_t, lspg_jacobian_policy_t>;

// };//end class



// template <typename fom_type,
// 	  ode::ImplicitEnum odeName,
// 	  typename decoder_type,
// 	  typename lspg_state_type>
// struct MaskedLSPGTypeGenerator
//   : LSPGCommonTypes<fom_type, odeName, decoder_type, lspg_state_type>{

//   using base_t = LSPGCommonTypes<fom_type, odeName, decoder_type, lspg_state_type>;

//   using typename base_t::fom_t;
//   using typename base_t::scalar_t;
//   using typename base_t::fom_state_t;
//   using typename base_t::fom_state_w_t;
//   using typename base_t::fom_rhs_w_t;
//   using typename base_t::decoder_t;
//   using typename base_t::decoder_jac_t;
//   using typename base_t::lspg_state_t;
//   using typename base_t::lspg_residual_t;
//   using typename base_t::fom_states_data;
//   using typename base_t::fom_rhs_data;

//   /* lspg_matrix_t is type of J*decoder_jac_t (in the most basic case) where
//    * * J is the jacobian of the fom rhs
//    * * decoder_jac_t is the type of the decoder jacobian
//    * In more complex cases, we might have (something)*J*decoder_jac_t,
//    * where (something) is product of few matrices.
//    * For now, set lspg_matrix_t to be of same type as decoder_jac_t
//    * if phi is MV<>, then lspg_matrix_t = core::MV<>
//    * if phi is Matrix<>, then we have core::Matrix<>
//    * not a bad assumption since all matrices are left-applied to decoder_jac_t
//    */
//   using lspg_matrix_t		= decoder_jac_t;

//   // policy for evaluating the rhs of the fom object
//   using fom_eval_rhs_policy_t	= rom::policy::EvaluateFomRhsDefault;

//   // policy for left multiplying the fom jacobian with decoder_jac_t
//   // possibly involving other stuff like explained above
//   using fom_apply_jac_policy_t	= rom::policy::ApplyFomJacobianDefault;

//   // policy defining how to compute the LSPG time-discrete residual
//   using lspg_residual_policy_t =
//     rom::decorator::Masked<
//     rom::LSPGResidualPolicy<
//       fom_states_data, fom_rhs_data, fom_eval_rhs_policy_t
//       >
//     >;

//   // policy defining how to compute the LSPG time-discrete jacobian
//   using lspg_jacobian_policy_t	=
//     rom::decorator::Masked<
//     rom::LSPGJacobianPolicy<
//       fom_states_data, lspg_matrix_t, fom_apply_jac_policy_t
//       >
//     >;

//   using aux_stepper_t = typename auxStepperHelper<odeName, lspg_state_type,
// 						  lspg_residual_t, lspg_matrix_t,
// 						  fom_type, lspg_residual_policy_t,
// 						  lspg_jacobian_policy_t>::type;

//   // declare type of stepper object
//   using rom_stepper_t		= ode::ImplicitStepper<
//     odeName, lspg_state_type, lspg_residual_t, lspg_matrix_t,
//     fom_type, aux_stepper_t, lspg_residual_policy_t, lspg_jacobian_policy_t>;

// };//end class


}}//end  namespace rompp::rom
#endif
