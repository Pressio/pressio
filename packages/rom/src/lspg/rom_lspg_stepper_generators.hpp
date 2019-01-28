
#ifndef ROM_LSPG_STEPPER_GENERATORS_HPP_
#define ROM_LSPG_STEPPER_GENERATORS_HPP_

#include "rom_lspg_type_generators.hpp"

namespace rompp{ namespace rom{

template <typename lspg_type_generator_t>
struct StepperObjectGenerator;


template <typename fom_type,
	  ode::ImplicitEnum odeName,
	  typename decoder_type,
	  typename lspg_state_type>
struct StepperObjectGenerator<
  PreconditionedLSPGTypeGenerator<fom_type, odeName, decoder_type, lspg_state_type>
  >{

  using types_info = PreconditionedLSPGTypeGenerator<
			fom_type, odeName, decoder_type, lspg_state_type>;
  using lspg_common		= typename types_info::base_t;

  using fom_t			= typename lspg_common::fom_t;
  using scalar_t		= typename lspg_common::scalar_t;
  using lspg_state_t		= typename lspg_common::lspg_state_t;
  using fom_state_t		= typename lspg_common::fom_state_t;
  using fom_state_w_t		= typename lspg_common::fom_state_w_t;
  using fom_rhs_w_t		= typename lspg_common::fom_rhs_w_t;
  using decoder_t		= typename lspg_common::decoder_t;
  using fom_states_data		= typename lspg_common::fom_states_data;
  using fom_rhs_data		= typename lspg_common::fom_rhs_data;

  using lspg_matrix_t		= typename types_info::lspg_matrix_t;
  using lspg_residual_policy_t	= typename types_info::lspg_residual_policy_t;
  using lspg_jacobian_policy_t	= typename types_info::lspg_jacobian_policy_t;
  using fom_eval_rhs_policy_t	= typename types_info::fom_eval_rhs_policy_t;
  using fom_apply_jac_policy_t	= typename types_info::fom_apply_jac_policy_t;
  using rom_stepper_t		= typename types_info::rom_stepper_t;

  fom_state_w_t y0Fom_		 = {};
  fom_rhs_w_t r0Fom_		 = {};
  fom_states_data fomStates_	 = {};
  fom_rhs_data fomRhs_		 = {};

  lspg_matrix_t romMat_		 = {};
  lspg_residual_policy_t resObj_ = {};
  lspg_jacobian_policy_t jacObj_ = {};
  rom_stepper_t stepperObj_	 = {};

  StepperObjectGenerator(const fom_t	   & appObj,
			 const fom_state_t & y0n,
			 decoder_t	   & decoder,
			 lspg_state_t	   & yROM,
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



template <typename fom_type,
	  ode::ImplicitEnum odeName,
	  typename decoder_type,
	  typename lspg_state_type>
struct StepperObjectGenerator<
  DefaultLSPGTypeGenerator<fom_type, odeName, decoder_type, lspg_state_type>
  >{

  using types_info = DefaultLSPGTypeGenerator<
			fom_type, odeName, decoder_type, lspg_state_type>;
  using lspg_common		= typename types_info::base_t;

  using fom_t			= typename lspg_common::fom_t;
  using scalar_t		= typename lspg_common::scalar_t;
  using lspg_state_t		= typename lspg_common::lspg_state_t;
  using fom_state_t		= typename lspg_common::fom_state_t;
  using fom_state_w_t		= typename lspg_common::fom_state_w_t;
  using fom_rhs_w_t		= typename lspg_common::fom_rhs_w_t;
  using decoder_t		= typename lspg_common::decoder_t;
  using fom_states_data		= typename lspg_common::fom_states_data;
  using fom_rhs_data		= typename lspg_common::fom_rhs_data;

  using lspg_matrix_t		= typename types_info::lspg_matrix_t;
  using lspg_residual_policy_t	= typename types_info::lspg_residual_policy_t;
  using lspg_jacobian_policy_t	= typename types_info::lspg_jacobian_policy_t;
  using fom_eval_rhs_policy_t	= typename types_info::fom_eval_rhs_policy_t;
  using fom_apply_jac_policy_t	= typename types_info::fom_apply_jac_policy_t;
  using rom_stepper_t		= typename types_info::rom_stepper_t;

  fom_state_w_t y0Fom_		 = {};
  fom_rhs_w_t r0Fom_		 = {};
  fom_states_data fomStates_	 = {};
  fom_rhs_data fomRhs_		 = {};

  lspg_matrix_t romMat_		 = {};
  lspg_residual_policy_t resObj_ = {};
  lspg_jacobian_policy_t jacObj_ = {};
  rom_stepper_t stepperObj_	 = {};

  StepperObjectGenerator(const fom_t	   & appObj,
			 const fom_state_t & y0n,
			 decoder_t	   & decoder,
			 lspg_state_t	   & yROM,
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
