
#ifndef ROM_GALERKIN_STEPPER_GENERATORS_HPP_
#define ROM_GALERKIN_STEPPER_GENERATORS_HPP_

#include "rom_galerkin_type_generators.hpp"

namespace rompp{ namespace rom{

template <typename T>
struct StepperObjectGenerator<
  T, core::meta::enable_if_t<
       std::is_void<typename T::aux_stepper_t>::value
       >
  > : T {

  using typename T::base_t::fom_t;
  using typename T::base_t::scalar_t;
  using typename T::base_t::rom_state_t;
  using typename T::base_t::fom_state_t;
  using typename T::base_t::fom_state_w_t;
  using typename T::base_t::fom_rhs_w_t;
  using typename T::base_t::decoder_t;
  using typename T::base_t::fom_states_data;
  using typename T::base_t::fom_rhs_data;

  using typename T::rom_matrix_t;
  using typename T::fom_eval_rhs_policy_t;
  using typename T::fom_apply_jac_policy_t;
  using typename T::rom_residual_policy_t;
  using typename T::rom_jacobian_policy_t;
  using typename T::rom_stepper_t;
  using typename T::aux_stepper_t;

  fom_eval_rhs_policy_t rhsEv_;   //default constructed
  fom_apply_jac_policy_t ajacEv_; //default constructed

  fom_state_w_t y0Fom_		 = {};
  fom_rhs_w_t r0Fom_		 = {};
  fom_states_data fomStates_	 = {};
  fom_rhs_data fomRhs_		 = {};

  rom_matrix_t romMat_		 = {};
  rom_residual_policy_t resObj_ = {};
  rom_jacobian_policy_t jacObj_ = {};
  rom_stepper_t stepperObj_	 = {};

  StepperObjectGenerator(const fom_t	   & appObj,
			 const fom_state_t & y0n,
			 decoder_t	   & decoder,
			 rom_state_t	   & yROM,
			 scalar_t	   t0)
    : y0Fom_(y0n),
      r0Fom_(rhsEv_.evaluate(appObj, y0Fom_, t0)),
      fomStates_(y0Fom_, decoder),
      fomRhs_(r0Fom_),
      romMat_(ajacEv_.evaluate(appObj, y0Fom_,
			       decoder.getJacobianRef(), t0)),
      resObj_(fomStates_, fomRhs_, rhsEv_),
      jacObj_(fomStates_, ajacEv_, romMat_),
      stepperObj_(appObj, resObj_, jacObj_, yROM)
  {}
};



template <typename T>
struct StepperObjectGenerator<
  T, core::meta::enable_if_t<
       !std::is_void<typename T::aux_stepper_t>::value
       >
  > : T {

  using typename T::base_t::fom_t;
  using typename T::base_t::scalar_t;
  using typename T::base_t::rom_state_t;
  using typename T::base_t::fom_state_t;
  using typename T::base_t::fom_state_w_t;
  using typename T::base_t::fom_rhs_w_t;
  using typename T::base_t::decoder_t;
  using typename T::base_t::fom_states_data;
  using typename T::base_t::fom_rhs_data;

  using typename T::rom_matrix_t;
  using typename T::fom_eval_rhs_policy_t;
  using typename T::fom_apply_jac_policy_t;
  using typename T::rom_residual_policy_t;
  using typename T::rom_jacobian_policy_t;
  using typename T::rom_stepper_t;
  using typename T::aux_stepper_t;

  fom_state_w_t y0Fom_		 = {};
  fom_rhs_w_t r0Fom_		 = {};
  fom_states_data fomStates_	 = {};
  fom_rhs_data fomRhs_		 = {};

  rom_matrix_t romMat_		 = {};
  rom_residual_policy_t resObj_ = {};
  rom_jacobian_policy_t jacObj_ = {};
  aux_stepper_t auxStepperObj_   = {};
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
      auxStepperObj_(appObj, resObj_, jacObj_, yROM),
      stepperObj_(appObj, resObj_, jacObj_, yROM, auxStepperObj_)
  {}
};

}}//end namespace rompp::rom
#endif
