
#ifndef ROM_DEFAULT_GALERKIN_EXPLICIT_VELOCITY_POLICY_HPP_
#define ROM_DEFAULT_GALERKIN_EXPLICIT_VELOCITY_POLICY_HPP_

#include "../rom_fwd.hpp"
#include "../../../ode/src/explicit/policies/ode_explicit_velocity_policy_base.hpp"
#include "../rom_data_fom_states.hpp"

namespace pressio{ namespace rom{

template <
  typename fom_states_data,
  typename fom_velocity_data,
  typename decoder_t
  >
class DefaultGalerkinExplicitVelocityPolicy
  : public ode::policy::ExplicitVelocityPolicyBase<
       DefaultGalerkinExplicitVelocityPolicy<fom_states_data,
					     fom_velocity_data,
					     decoder_t>>,
    protected fom_states_data,
    protected fom_velocity_data{

protected:
  using this_t = DefaultGalerkinExplicitVelocityPolicy<
		      fom_states_data, fom_velocity_data, decoder_t>;
  using base_t = ode::policy::ExplicitVelocityPolicyBase<this_t>;
  friend base_t;

  const decoder_t & decoder_;
  using fom_states_data::yFom_;
  using fom_velocity_data::fomRhs_;

public:
  static constexpr bool isResidualPolicy_ = true;

public:
  DefaultGalerkinExplicitVelocityPolicy() = delete;
  ~DefaultGalerkinExplicitVelocityPolicy() = default;
  DefaultGalerkinExplicitVelocityPolicy(const fom_states_data & fomStates,
					const fom_velocity_data & fomResids,
					const decoder_t & decoder)
    : fom_states_data(fomStates),
      fom_velocity_data(fomResids),
      decoder_(decoder){}

public:
  /* for now, the state and velocity must be of the same type */

  template <
    typename galerkin_state_t,
    typename fom_t,
    typename scalar_t
  >
  void operator()(const galerkin_state_t  & romY,
		  galerkin_state_t	  & romR,
  		  const fom_t		  & app,
		  scalar_t		  t) const
  {
#ifdef HAVE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("galerkin explicit velocity");
#endif

    fom_states_data::template reconstructCurrentFomState(romY);

#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("fom eval rhs");
#endif
    app.velocity(*yFom_.data(), *fomRhs_.data(), t);
#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("fom eval rhs");
#endif

#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("phiT*fomRhs");
#endif
    containers::ops::dot(decoder_.getReferenceToJacobian(), fomRhs_, romR);
#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("phiT*fomRhs");
#endif

#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("galerkin explicit velocity");
#endif
  }


  template <
    typename galerkin_state_t,
    typename fom_t,
    typename scalar_t
    >
  galerkin_state_t operator()(const galerkin_state_t  & romY,
			      const fom_t		 & app,
			      scalar_t		 t) const
  {
    // for now, make it better later
    galerkin_state_t result(romY);
    result.setZero();

#ifdef HAVE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("galerkin explicit velocity");
#endif

    fom_states_data::template reconstructCurrentFomState(romY);

#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("fom eval rhs");
#endif
    typename fom_velocity_data::fom_rhs_t fomR(app.velocity(*yFom_.data(), t));
#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("fom eval rhs");
#endif

#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("phiT*fomRhs");
#endif
    containers::ops::dot(decoder_.getReferenceToJacobian(), fomR, result);
#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("phiT*fomRhs");
#endif

#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("galerkin explicit velocity");
#endif

    return result;
  }

};//end class

}}//end namespace pressio::rom
#endif
