
#ifndef ROM_LSPG_STEPPER_GENERATORS_HPP_
#define ROM_LSPG_STEPPER_GENERATORS_HPP_

#include "rom_lspg_type_generators.hpp"

namespace rompp{ namespace rom{

template <typename lspg_problem>
struct LSPGStepperObjectGenerator<
  lspg_problem,
  core::meta::enable_if_t<
    std::is_void<typename lspg_problem::aux_stepper_t>::value
    >
  > : lspg_problem {

  using typename lspg_problem::fom_t;
  using typename lspg_problem::scalar_t;
  using typename lspg_problem::lspg_state_t;
  using typename lspg_problem::fom_state_t;
  using typename lspg_problem::fom_state_w_t;
  using typename lspg_problem::fom_rhs_w_t;
  using typename lspg_problem::decoder_t;
  using typename lspg_problem::fom_states_data;
  using typename lspg_problem::fom_rhs_data;

  using typename lspg_problem::lspg_matrix_t;
  using typename lspg_problem::fom_eval_rhs_policy_t;
  using typename lspg_problem::fom_apply_jac_policy_t;
  using typename lspg_problem::lspg_residual_policy_t;
  using typename lspg_problem::lspg_jacobian_policy_t;
  using typename lspg_problem::rom_stepper_t;
  using typename lspg_problem::aux_stepper_t;

  fom_eval_rhs_policy_t rhsEv_;
  fom_apply_jac_policy_t ajacEv_;

  fom_state_w_t y0Fom_;
  fom_rhs_w_t r0Fom_;
  fom_states_data fomStates_;
  fom_rhs_data fomRhs_;
  lspg_matrix_t romMat_;
  lspg_residual_policy_t resObj_;
  lspg_jacobian_policy_t jacObj_;
  rom_stepper_t stepperObj_;

  LSPGStepperObjectGenerator(const fom_t	& appObj,
			     const fom_state_t	& y0n,
			     decoder_t		& decoder,
			     lspg_state_t	& yROM,
			     scalar_t		t0)
    : y0Fom_(y0n),
      r0Fom_(rhsEv_.evaluate(appObj, y0Fom_, t0)),
      fomStates_(y0Fom_, decoder),
      fomRhs_(r0Fom_),
      romMat_(ajacEv_.evaluate(appObj, y0Fom_,
			       decoder.getJacobianRef(), t0)),
      resObj_(fomStates_, fomRhs_, rhsEv_),
      jacObj_(fomStates_, ajacEv_, romMat_),
      stepperObj_(appObj, resObj_, jacObj_, yROM){}

};



template <typename lspg_problem>
struct LSPGStepperObjectGenerator<
  lspg_problem, core::meta::enable_if_t<
       !std::is_void<typename lspg_problem::aux_stepper_t>::value
       >
  > : lspg_problem {

  using typename lspg_problem::fom_t;
  using typename lspg_problem::scalar_t;
  using typename lspg_problem::lspg_state_t;
  using typename lspg_problem::fom_state_t;
  using typename lspg_problem::fom_state_w_t;
  using typename lspg_problem::fom_rhs_w_t;
  using typename lspg_problem::decoder_t;
  using typename lspg_problem::fom_states_data;
  using typename lspg_problem::fom_rhs_data;

  using typename lspg_problem::lspg_matrix_t;
  using typename lspg_problem::fom_eval_rhs_policy_t;
  using typename lspg_problem::fom_apply_jac_policy_t;
  using typename lspg_problem::lspg_residual_policy_t;
  using typename lspg_problem::lspg_jacobian_policy_t;
  using typename lspg_problem::rom_stepper_t;
  using typename lspg_problem::aux_stepper_t;

  fom_eval_rhs_policy_t rhsEv_;
  fom_apply_jac_policy_t ajacEv_;

  fom_state_w_t y0Fom_;
  fom_rhs_w_t r0Fom_;
  fom_states_data fomStates_;
  fom_rhs_data fomRhs_;
  lspg_matrix_t romMat_;
  lspg_residual_policy_t resObj_;
  lspg_jacobian_policy_t jacObj_;
  aux_stepper_t auxStepperObj_;
  rom_stepper_t stepperObj_;

  LSPGStepperObjectGenerator(const fom_t	& appObj,
			     const fom_state_t  & y0n,
			     decoder_t		& decoder,
			     lspg_state_t	& yROM,
			     scalar_t		t0)
    : y0Fom_(y0n),
      r0Fom_(rhsEv_.evaluate(appObj, y0Fom_, t0)),
      fomStates_(y0Fom_, decoder),
      fomRhs_(r0Fom_),
      romMat_(ajacEv_.evaluate(appObj, y0Fom_,
			       decoder.getJacobianRef(), t0)),
      resObj_(fomStates_, fomRhs_, rhsEv_),
      jacObj_(fomStates_, ajacEv_, romMat_),
      auxStepperObj_(appObj, resObj_, jacObj_, yROM),
      stepperObj_(appObj, resObj_, jacObj_, yROM, auxStepperObj_)
  {}
};

}}//end namespace rompp::rom
#endif
