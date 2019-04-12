
#ifndef ROM_GALERKIN_STEPPER_GENERATORS_HPP_
#define ROM_GALERKIN_STEPPER_GENERATORS_HPP_

#include "rom_galerkin_type_generators.hpp"

namespace rompp{ namespace rom{

template <typename problem_types>
struct GalerkinStepperObjectGenerator<
  problem_types,
  ::rompp::mpl::enable_if_t<
    problem_types::odeName_ == ode::ExplicitEnum::Euler or
    problem_types::odeName_ == ode::ExplicitEnum::RungeKutta4
    >
  > : problem_types {

  using typename problem_types::base_t::fom_t;
  using typename problem_types::base_t::scalar_t;
  using typename problem_types::base_t::galerkin_state_t;
  using typename problem_types::base_t::galerkin_residual_t;
  using typename problem_types::base_t::fom_state_t;
  using typename problem_types::base_t::fom_state_w_t;
  using typename problem_types::base_t::fom_rhs_w_t;
  using typename problem_types::base_t::decoder_t;
  using typename problem_types::base_t::fom_states_data;
  using typename problem_types::base_t::fom_rhs_data;

  using typename problem_types::galerkin_residual_policy_t;
  using typename problem_types::galerkin_stepper_t;

  fom_state_w_t y0Fom_			= {};
  fom_rhs_w_t r0Fom_			= {};
  fom_states_data fomStates_		= {};
  fom_rhs_data fomRhs_			= {};
  galerkin_residual_policy_t resObj_	= {};
  galerkin_stepper_t stepperObj_	= {};

  GalerkinStepperObjectGenerator(const fom_t	   & appObj,
				 const fom_state_t & y0n,
				 decoder_t	   & decoder,
				 galerkin_state_t  & yROM,
				 scalar_t	   t0)
    : y0Fom_(y0n),
      r0Fom_(appObj.residual(y0n, t0)),
      fomStates_(y0Fom_, decoder),
      fomRhs_(r0Fom_),
      resObj_(fomStates_, fomRhs_, decoder),
      stepperObj_(appObj, resObj_, yROM)
  {}

};

}}//end namespace rompp::rom
#endif
