
#ifndef ROM_GALERKIN_STEPPER_GENERATORS_HPP_
#define ROM_GALERKIN_STEPPER_GENERATORS_HPP_

#include "rom_galerkin_type_generators.hpp"

namespace rompp{ namespace rom{

template <typename T>
struct GalerkinExplicitStepperObjectGenerator<T> : T {

  using typename T::base_t::fom_t;
  using typename T::base_t::scalar_t;
  using typename T::base_t::galerkin_state_t;
  using typename T::base_t::galerkin_residual_t;
  using typename T::base_t::fom_state_t;
  using typename T::base_t::fom_state_w_t;
  using typename T::base_t::fom_rhs_w_t;
  using typename T::base_t::decoder_t;
  using typename T::base_t::fom_states_data;
  using typename T::base_t::fom_rhs_data;

  using typename T::galerkin_residual_policy_t;
  using typename T::galerkin_stepper_t;

  fom_state_w_t y0Fom_		 = {};
  fom_rhs_w_t r0Fom_		 = {};
  fom_states_data fomStates_	 = {};
  fom_rhs_data fomRhs_		 = {};

  galerkin_residual_policy_t resObj_ = {};
  galerkin_stepper_t stepperObj_	 = {};

  GalerkinExplicitStepperObjectGenerator(const fom_t	   & appObj,
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
