
#ifndef ROM_GALERKIN_STEPPER_GENERATORS_HPP_
#define ROM_GALERKIN_STEPPER_GENERATORS_HPP_

#include "rom_galerkin_type_generators.hpp"

namespace pressio{ namespace rom{

template <typename problem_t>
struct GalerkinProblemGenerator<problem_t>
  : problem_t {

  using typename problem_t::fom_t;
  using typename problem_t::scalar_t;
  using typename problem_t::fom_native_state_t;
  using typename problem_t::fom_state_t;
  using typename problem_t::fom_velocity_t;

  using typename problem_t::galerkin_state_t;
  using typename problem_t::decoder_t;
  using typename problem_t::fom_state_reconstr_t;
  using typename problem_t::fom_states_data;
  using typename problem_t::fom_velocity_data;

  using typename problem_t::galerkin_residual_policy_t;
  using typename problem_t::galerkin_stepper_t;

  fom_state_t			yFomRef_;
  fom_state_reconstr_t		yFomReconstructor_;
  fom_velocity_t			rFomRef_;
  fom_states_data		fomStates_;
  fom_velocity_data			fomRhs_;
  galerkin_residual_policy_t	resObj_;
  galerkin_stepper_t		stepperObj_;

  GalerkinProblemGenerator(const fom_t		    & appObj,
  			   const fom_native_state_t & yFomRefNative,
  			   decoder_t		    & decoder,
  			   galerkin_state_t	    & yROM,
  			   scalar_t		    t0)
    : yFomRef_(yFomRefNative),
      yFomReconstructor_(yFomRef_, decoder),
      rFomRef_( appObj.velocity(*yFomRef_.data(), t0) ),
      fomStates_(yFomRef_, yFomReconstructor_),
      fomRhs_(rFomRef_),
      resObj_(fomStates_, fomRhs_, decoder),
      stepperObj_(yROM, appObj, resObj_)
  {}

};

}}//end namespace pressio::rom
#endif
