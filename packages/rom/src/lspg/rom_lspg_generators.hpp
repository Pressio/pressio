
#ifndef ROM_LSPG_TYPE_GENERATOR_HPP_
#define ROM_LSPG_TYPE_GENERATOR_HPP_

#include "../rom_ConfigDefs.hpp"
#include "../rom_forward_declarations.hpp"
#include "../rom_data_fom_rhs.hpp"
#include "../rom_data_fom_states.hpp"
#include "../policies/rom_evaluate_fom_rhs_policy.hpp"
#include "../policies/rom_apply_fom_jacobian_policy.hpp"

namespace rompp{ namespace rom{

template <ode::ImplicitEnum odeName>
struct statesStorageHelper;

template <>
struct statesStorageHelper<ode::ImplicitEnum::Euler>{
  static constexpr int maxAuxStates_ = 1;
};

template <>
struct statesStorageHelper<ode::ImplicitEnum::BDF2>{
  static constexpr int maxAuxStates_ = 2;
};
//-------------------------------------------------------


template <typename fom_type,
	  ode::ImplicitEnum odeName,
	  typename decoder_type,
	  typename rom_state_type>
struct DefaultLSPGTypeGenerator{

  // these are native types of the full-order model (fom)
  using fom_t			= fom_type;
  using scalar_t		= typename fom_t::scalar_type;
  using fom_state_t		= typename fom_t::state_type;
  using fom_rhs_t		= typename fom_t::residual_type;

  // declare fom wrapper types
  using fom_state_w_t		= rompp::core::Vector<fom_state_t>;
  using fom_rhs_w_t		= rompp::core::Vector<fom_rhs_t>;

  // decoder types (passed in)
  using decoder_t		= decoder_type;
  using decoder_jac_t		= typename decoder_t::jacobian_t;

  // rom state type (passed in)
  using rom_state_t		= rom_state_type;

  // for LSPG, the rom residual type = core::wrapper of application rhs
  // i.e. the wrapped fom rhs type
  using rom_residual_t		= fom_rhs_w_t;

  // max num of states needed for time integration.
  // this is deduced based on odeName
  static constexpr auto maxAuxStates = statesStorageHelper<odeName>::maxAuxStates_;

  // class type holding fom states data
  using fom_states_data = rompp::rom::FomStatesData<
	fom_state_w_t, maxAuxStates, decoder_t>;

  // class type holding fom rhs data
  using fom_rhs_data = rompp::rom::FomRhsData<fom_rhs_w_t>;

  /* rom_matrix_t is type of J*decoder_jac_t (in the most basic case) where
   * * J is the jacobian of the fom rhs
   * * decoder_jac_t is the type of the decoder jacobian
   * In more complex cases, we might have (something)*J*decoder_jac_t,
   * where (something) is product of few matrices.
   * For now, set rom_matrix_t to be of same type as decoder_jac_t
   * if phi is MV<>, then rom_matrix_t = core::MV<>
   * if phi is Matrix<>, then we have core::Matrix<>
   * not a bad assumption since all matrices are left-applied to decoder_jac_t
   */
  using rom_matrix_t		= decoder_jac_t;

  // policy for evaluating the rhs of the fom object
  using fom_eval_rhs_policy_t	= rompp::rom::policy::EvaluateFomRhsDefault;

  // policy for left multiplying the fom jacobian with decoder_jac_t
  // possibly involving other stuff like explained above
  using fom_apply_jac_policy_t	= rompp::rom::policy::ApplyFomJacobianDefault;

  // policy defining how to compute the LSPG time-discrete residual
  using lspg_residual_policy_t	= rompp::rom::LSPGResidualPolicy<
	fom_states_data, fom_rhs_data, fom_eval_rhs_policy_t>;

  // policy defining how to compute the LSPG time-discrete jacobian
  using lspg_jacobian_policy_t	= rompp::rom::LSPGJacobianPolicy<
	fom_states_data, rom_matrix_t, fom_apply_jac_policy_t>;

  // declare type of stepper object
  using rom_stepper_t		= rompp::ode::ImplicitStepper<
    odeName, rom_state_t, rom_residual_t, rom_matrix_t,
    fom_t, void, lspg_residual_policy_t, lspg_jacobian_policy_t>;

};//end class




template <typename lspg_type_info>
struct StepperObjectGenerator;


template <typename fom_type,
	  ode::ImplicitEnum odeName,
	  typename decoder_type,
	  typename rom_state_type>
struct StepperObjectGenerator<
  DefaultLSPGTypeGenerator<fom_type, odeName, decoder_type, rom_state_type>
  >{
  using types_info = DefaultLSPGTypeGenerator<
    fom_type, odeName, decoder_type, rom_state_type>;

  using fom_t			= typename types_info::fom_t;
  using scalar_t		= typename types_info::scalar_t;
  using rom_state_t		= typename types_info::rom_state_t;
  using fom_state_t		= typename types_info::fom_state_t;
  using fom_state_w_t		= typename types_info::fom_state_w_t;
  using fom_rhs_w_t		= typename types_info::fom_rhs_w_t;
  using decoder_t		= typename types_info::decoder_t;
  using rom_matrix_t		= typename types_info::rom_matrix_t;
  using fom_states_data		= typename types_info::fom_states_data;
  using fom_rhs_data		= typename types_info::fom_rhs_data;
  using lspg_residual_policy_t	= typename types_info::lspg_residual_policy_t;
  using lspg_jacobian_policy_t	= typename types_info::lspg_jacobian_policy_t;
  using fom_eval_rhs_policy_t	= typename types_info::fom_eval_rhs_policy_t;
  using fom_apply_jac_policy_t	= typename types_info::fom_apply_jac_policy_t;
  using rom_stepper_t		= typename types_info::rom_stepper_t;

  fom_state_w_t y0Fom_		 = {};
  fom_rhs_w_t r0Fom_		 = {};
  fom_states_data fomStates_	 = {};
  fom_rhs_data fomRhs_		 = {};

  rom_matrix_t romMat_		 = {};
  lspg_residual_policy_t resObj_ = {};
  lspg_jacobian_policy_t jacObj_ = {};
  rom_stepper_t stepperObj_	 = {};

  StepperObjectGenerator(const fom_t	   & appObj,
			 const fom_state_t & y0n,
			 decoder_t	   & decoder,
			 rom_state_t	   & yROM,
			 scalar_t	   t0)
    : y0Fom_(y0n),
      r0Fom_(fom_eval_rhs_policy_t().evaluate(appObj, y0Fom_, t0)),
      fomStates_(y0Fom_, decoder),
      fomRhs_(r0Fom_),
      romMat_(fom_apply_jac_policy_t().evaluate(appObj, y0Fom_,
						decoder.getJacobianRef(),
						t0)),
      resObj_(fomStates_, fomRhs_, fom_eval_rhs_policy_t()),
      jacObj_(fomStates_, fom_apply_jac_policy_t(), romMat_),
      stepperObj_(appObj, resObj_, jacObj_, yROM)
  {}
};


}}//end namespace rompp::rom
#endif
