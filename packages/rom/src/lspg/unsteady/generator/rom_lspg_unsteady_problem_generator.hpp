
#ifndef ROMPP_ROM_LSPG_UNSTEADY_PROBLEM_GENERATOR_HPP_
#define ROMPP_ROM_LSPG_UNSTEADY_PROBLEM_GENERATOR_HPP_

#include "rom_lspg_type_generator_default.hpp"
#include "rom_lspg_type_generator_preconditioned.hpp"
#include "rom_lspg_type_generator_masked.hpp"

namespace rompp{ namespace rom{

template <typename lspg_problem>
struct LSPGUnsteadyProblemGenerator<
  lspg_problem> : lspg_problem {

  using typename lspg_problem::fom_t;
  using typename lspg_problem::scalar_t;
  using typename lspg_problem::fom_state_t;
  using typename lspg_problem::fom_state_w_t;
  using typename lspg_problem::fom_rhs_w_t;

  using typename lspg_problem::lspg_state_t;
  using typename lspg_problem::decoder_t;
  using typename lspg_problem::fom_state_reconstr_t;
  using typename lspg_problem::fom_states_data;
  using typename lspg_problem::fom_rhs_data;
  using typename lspg_problem::td_ud_ops;

  using typename lspg_problem::lspg_matrix_t;
  using typename lspg_problem::fom_eval_rhs_policy_t;
  using typename lspg_problem::fom_apply_jac_policy_t;
  using typename lspg_problem::lspg_residual_policy_t;
  using typename lspg_problem::lspg_jacobian_policy_t;

  using typename lspg_problem::aux_stepper_t;
  using typename lspg_problem::lspg_stepper_t;

  fom_eval_rhs_policy_t		rhsEv_;
  fom_apply_jac_policy_t	ajacEv_;
  fom_state_w_t			yFomRef_;
  fom_state_reconstr_t		yFomReconstructor_;
  fom_rhs_w_t			rFomRef_;
  fom_states_data		fomStates_;
  fom_rhs_data			fomRhs_;
  lspg_matrix_t			romMat_;
  const td_ud_ops	      & tdOps_;
  lspg_residual_policy_t	resObj_;
  lspg_jacobian_policy_t	jacObj_;

  /* here we use conditional type if auxiliary stepper is non-void,
   * otherwise we set it to a dummy type and we dont construct it */
  typename std::conditional<
    std::is_void<aux_stepper_t>::value,
    core::impl::empty, aux_stepper_t>::type auxStepperObj_;

  // actual stepper object
  lspg_stepper_t			stepperObj_;

public:
  /* sfinae here for when we do NOT need aux stepper
   * note that we need to use trick _fom_t for sfinae to work */
  template <
   typename _fom_t,
   typename ::rompp::mpl::enable_if_t<
     std::is_void<aux_stepper_t>::value and
     std::is_same<_fom_t, fom_t>::value
     > * = nullptr
  >
  LSPGUnsteadyProblemGenerator(const _fom_t	 & appObj,
			       const fom_state_t & yFomRefNative,
			       decoder_t	 & decoder,
			       lspg_state_t	 & yROM,
			       scalar_t		 t0,
			       const td_ud_ops & tdOps)
    : rhsEv_{},
      ajacEv_{},
      yFomRef_(yFomRefNative),
      yFomReconstructor_(yFomRef_, decoder),
      rFomRef_( rhsEv_.evaluate(appObj, yFomRef_, t0) ),
      fomStates_(yFomRef_, yFomReconstructor_),
      fomRhs_(rFomRef_),
      romMat_(ajacEv_.evaluate(appObj, yFomRef_,
			       decoder.getReferenceToJacobian(), t0)),
      tdOps_{tdOps},
      resObj_(fomStates_, fomRhs_, rhsEv_, tdOps_),
      jacObj_(fomStates_, ajacEv_, romMat_, decoder, tdOps_),
      auxStepperObj_{},
      stepperObj_(yROM, appObj, resObj_, jacObj_)
  {}

  /* sfinae here for when we need aux stepper
   * note that we need to use trick _fom_t for sfinae to work */
  template <
    typename _fom_t,
    typename ::rompp::mpl::enable_if_t<
      std::is_void<aux_stepper_t>::value == false and
      std::is_same<_fom_t, fom_t>::value
      > * = nullptr
    >
  LSPGUnsteadyProblemGenerator(const _fom_t	 & appObj,
			       const fom_state_t & yFomRefNative,
			       const decoder_t	 & decoder,
			       lspg_state_t	 & yROM,
			       scalar_t		 t0)
    : rhsEv_{},
      ajacEv_{},
      yFomRef_(yFomRefNative),
      yFomReconstructor_(yFomRef_, decoder),
      rFomRef_( rhsEv_.evaluate(appObj, yFomRef_, t0) ),
      fomStates_(yFomRef_, yFomReconstructor_),
      fomRhs_(rFomRef_),
      romMat_(ajacEv_.evaluate(appObj, yFomRef_,
  			       decoder.getReferenceToJacobian(), t0)),
      resObj_(fomStates_, fomRhs_, rhsEv_),
      jacObj_(fomStates_, ajacEv_, romMat_, decoder),
      auxStepperObj_(yROM, appObj, resObj_, jacObj_),
      stepperObj_(yROM, appObj, resObj_, jacObj_, auxStepperObj_)
  {}

};

}}//end namespace rompp::rom
#endif
