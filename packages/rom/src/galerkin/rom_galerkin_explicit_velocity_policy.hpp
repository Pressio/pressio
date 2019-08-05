
#ifndef ROM_DEFAULT_GALERKIN_EXPLICIT_VELOCITY_POLICY_HPP_
#define ROM_DEFAULT_GALERKIN_EXPLICIT_VELOCITY_POLICY_HPP_

#include "../rom_fwd.hpp"
#include "../../../ode/src/explicit/policies/ode_explicit_velocity_policy_base.hpp"
#include "../rom_data_fom_states.hpp"

namespace pressio{ namespace rom{

template <
  typename fom_states_data,
  typename fom_velocity_data,
  typename decoder_t,
  typename ud_ops
  >
class DefaultGalerkinExplicitVelocityPolicy
  : public ode::policy::ExplicitVelocityPolicyBase<
       DefaultGalerkinExplicitVelocityPolicy<fom_states_data,
					     fom_velocity_data,
					     decoder_t, ud_ops>>,
    protected fom_states_data,
    protected fom_velocity_data{

protected:
  using this_t = DefaultGalerkinExplicitVelocityPolicy<fom_states_data,
						       fom_velocity_data,
						       decoder_t,
						       ud_ops>;
  friend ode::policy::ExplicitVelocityPolicyBase<this_t>;

  const decoder_t & decoder_;
  using fom_states_data::yFom_;
  using fom_velocity_data::fomRhs_;

#ifdef HAVE_PYBIND11
  typename std::conditional<
    mpl::is_same<ud_ops, pybind11::object>::value,
    ud_ops,
    const ud_ops *
    >::type udOps_ = {};
#else
    const ud_ops * udOps_ = {};
#endif

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

#ifdef HAVE_PYBIND11
  // this cnstr only enabled when udOps is non-void and python
  template <
    typename _ud_ops = ud_ops,
    mpl::enable_if_t<
      mpl::is_same<_ud_ops, pybind11::object>::value
      > * = nullptr
    >
  DefaultGalerkinExplicitVelocityPolicy(const fom_states_data & fomStates,
					const fom_velocity_data & fomResids,
					const decoder_t & decoder,
					const _ud_ops & udOps)
    : fom_states_data(fomStates),
      fom_velocity_data(fomResids),
      decoder_(decoder),
      udOps_{udOps}{}
#endif

public:
  /* for now, the ROM state and ROM velocity must be of the same type */
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
    this->compute_impl(romY, romR, app, t);
  }

  template <
    typename galerkin_state_t,
    typename fom_t,
    typename scalar_t
    >
  galerkin_state_t operator()(const galerkin_state_t  & romY,
			      const fom_t	      & app,
			      scalar_t		 t) const
  {
    // for now, make it better later
    galerkin_state_t result(romY);
    ::pressio::containers::ops::set_zero(result);
    this->compute_impl(romY, result, app, t);
    return result;
  }

private:

  template <
  typename galerkin_state_t,
  typename fom_t,
  typename scalar_t,
  typename _ud_ops = ud_ops,
  mpl::enable_if_t<
    std::is_void<_ud_ops>::value and
    ::pressio::containers::meta::is_vector_wrapper<galerkin_state_t>::value
#ifdef HAVE_PYBIND11
    and !::pressio::containers::meta::is_cstyle_array_pybind11<galerkin_state_t>::value
#endif
    > * = nullptr
  >
  void compute_impl(const galerkin_state_t  & romY,
		    galerkin_state_t	    & romR,
		    const fom_t		    & app,
		    scalar_t		    t) const
  {
#ifdef HAVE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("galerkin explicit velocity");
#endif

    fom_states_data::template reconstructCurrentFomState(romY);

#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("fom eval rhs");
#endif
    app.velocity(*yFom_.data(), t, *fomRhs_.data());
#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("fom eval rhs");
#endif

#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("phiT*fomRhs");
#endif
    const auto & phi = decoder_.getReferenceToJacobian();
    containers::ops::dot(phi, fomRhs_, romR);
#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("phiT*fomRhs");
#endif

#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("galerkin explicit velocity");
#endif
  }



#ifdef HAVE_PYBIND11
  template <
  typename galerkin_state_t,
  typename fom_t,
  typename scalar_t,
  typename _ud_ops = ud_ops,
  mpl::enable_if_t<
    ::pressio::mpl::is_same<_ud_ops, pybind11::object>::value and
    ::pressio::containers::meta::is_cstyle_array_pybind11<galerkin_state_t>::value
    > * = nullptr
  >
  void compute_impl(const galerkin_state_t  & romY,
		    galerkin_state_t	    & romR,
		    const fom_t		    & app,
		    scalar_t		    t) const
  {
#ifdef HAVE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("galerkin explicit velocity");
#endif

    fom_states_data::template reconstructCurrentFomState(romY);

#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("fom eval rhs");
#endif
    app.attr("velocity")(yFom_, t, fomRhs_);
#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("fom eval rhs");
#endif

#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("phiT*fomRhs");
#endif
    const auto & phi = decoder_.getReferenceToJacobian();
    auto constexpr transposePhi = true;
    udOps_.attr("multiply2")(phi, fomRhs_, romR, transposePhi);

#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("phiT*fomRhs");
#endif

#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("galerkin explicit velocity");
#endif
  }
#endif

};//end class

}}//end namespace pressio::rom
#endif
