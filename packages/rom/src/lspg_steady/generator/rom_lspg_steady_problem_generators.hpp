
#ifndef ROM_LSPG_STEADY_PROBLEM_STEPPER_GENERATOR_HPP_
#define ROM_LSPG_STEADY_PROBLEM_STEPPER_GENERATOR_HPP_

#include "rom_lspg_type_generator_default.hpp"

namespace rompp{ namespace rom{

template <typename lspg_problem>
struct LSPGSteadyProblemGenerator<
  lspg_problem > : lspg_problem {

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
  using typename lspg_problem::rom_system_t;

  fom_eval_rhs_policy_t rhsEv_;
  fom_apply_jac_policy_t ajacEv_;

  fom_state_w_t y0Fom_;
  fom_rhs_w_t r0Fom_;
  fom_states_data fomStates_;
  fom_rhs_data fomRhs_;
  lspg_matrix_t romMat_;
  lspg_residual_policy_t resObj_;
  lspg_jacobian_policy_t jacObj_;
  rom_system_t systemObj_;

  LSPGSteadyProblemGenerator(const fom_t	& appObj,
			     const fom_state_t	& y0n,
			     decoder_t		& decoder,
			     lspg_state_t	& yROM)
    : y0Fom_(y0n),
      r0Fom_(rhsEv_.evaluate(appObj, y0Fom_)),
      fomStates_(y0Fom_, decoder),
      fomRhs_(r0Fom_),
      romMat_(ajacEv_.evaluate(appObj, y0Fom_,
			       decoder.getJacobianRef())),
      resObj_(fomStates_, fomRhs_, rhsEv_),
      jacObj_(fomStates_, ajacEv_, romMat_),
      systemObj_(appObj, resObj_, jacObj_){}

};

}}//end namespace rompp::rom
#endif
