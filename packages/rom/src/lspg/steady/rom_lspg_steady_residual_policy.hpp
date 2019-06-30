
#ifndef ROM_LSPG_STEADY_RESIDUAL_POLICY_HPP_
#define ROM_LSPG_STEADY_RESIDUAL_POLICY_HPP_

#include "../../rom_fwd.hpp"
#include "../../rom_data_fom_rhs.hpp"
#include "../../rom_data_fom_states.hpp"

namespace rompp{ namespace rom{

template <
  typename fom_states_data,
  typename fom_rhs_data,
  typename fom_eval_rhs_policy
  >
class LSPGSteadyResidualPolicy
  : protected fom_states_data,
    protected fom_rhs_data,
    protected fom_eval_rhs_policy{

protected:
  // protected because we might have decorators of this class
  using this_t = LSPGSteadyResidualPolicy<fom_states_data,
					  fom_rhs_data,
					  fom_eval_rhs_policy>;

  using fom_states_data::yFom_;
  using fom_rhs_data::fomRhs_;

public:
  static constexpr bool isResidualPolicy_ = true;
  using typename fom_rhs_data::fom_rhs_t;

public:
  LSPGSteadyResidualPolicy() = delete;
  ~LSPGSteadyResidualPolicy() = default;

  LSPGSteadyResidualPolicy(const fom_states_data     & fomStates,
			   const fom_rhs_data	     & fomResids,
			   const fom_eval_rhs_policy & fomEvalRhsFunctor)
    : fom_states_data(fomStates),
      fom_rhs_data(fomResids),
      fom_eval_rhs_policy(fomEvalRhsFunctor){}

public:
  template <typename lspg_state_t,
	    typename lspg_residual_t,
	    typename fom_t>
  void operator()(const lspg_state_t	& romY,
		  lspg_residual_t	& romR,
  		  const fom_t		& app) const
  {
#ifdef HAVE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("lspg residual");
#endif

    fom_states_data::template reconstructCurrentFomState(romY);

#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("fom eval rhs");
    fom_eval_rhs_policy::evaluate(app, yFom_, romR);
    timer->stop("fom eval rhs");
#else
    fom_eval_rhs_policy::evaluate(app, yFom_, romR);
#endif

#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("lspg residual");
#endif
  }

  template <typename lspg_state_t, typename fom_t>
  fom_rhs_t operator()(const lspg_state_t & romY,
			 const fom_t	    & app) const
  {
    (*this).template operator()(romY, fomRhs_, app);
    return fomRhs_;
  }

};//end class

}}//end namespace rompp::rom
#endif
