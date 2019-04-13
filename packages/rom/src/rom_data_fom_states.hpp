
#ifndef ROM_DATA_FOM_STATES_HPP_
#define ROM_DATA_FOM_STATES_HPP_

#include "rom_ConfigDefs.hpp"

namespace rompp{ namespace rom{

template <typename fom_state_w_type,
	  int maxNstates,
	  typename reconstuctor_type>
struct FomStatesData{

  static constexpr int maxNstates_ = maxNstates;
  using fom_state_w_t = fom_state_w_type;

  FomStatesData() = delete;
  ~FomStatesData() = default;

  template <
    int _maxNstates = maxNstates,
    ::rompp::mpl::enable_if_t<_maxNstates==0> * = nullptr
    >
  FomStatesData(const fom_state_w_t & yFomIn,
		const reconstuctor_type & fomStateReconstr)
    : yFom_(yFomIn),
      fomStateReconstrObj_(fomStateReconstr)
  {
    this->resetContainersToZero();
  }

  template <
    int _maxNstates = maxNstates,
    ::rompp::mpl::enable_if_t<_maxNstates==1> * = nullptr
    >
  FomStatesData(const fom_state_w_t & yFomIn,
		const reconstuctor_type & fomStateReconstr)
    : yFom_(yFomIn),
      yFomOld_{yFomIn},
      fomStateReconstrObj_(fomStateReconstr)
  {
    this->resetContainersToZero();
  }

  template <
    int _maxNstates = maxNstates,
    ::rompp::mpl::enable_if_t<_maxNstates==2> * = nullptr
    >
  FomStatesData(const fom_state_w_t & yFomIn,
		const reconstuctor_type & fomStateReconstr)
    : yFom_(yFomIn),
      yFomOld_{yFomIn, yFomIn},
      fomStateReconstrObj_(fomStateReconstr)
  {
    this->resetContainersToZero();
  }

protected:

  template <typename rom_state_t>
  void reconstructCurrentFomState(const rom_state_t & romY) const
  {
#ifdef HAVE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("reconstruct fom state");
#endif

    fomStateReconstrObj_(romY, yFom_);

#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("reconstruct fom state");
#endif
  }

  template <int n, typename rom_state_t>
  void reconstructFomOldStates(const std::array<
			       rom_state_t, n
			       > & romYprev) const
  {
#ifdef HAVE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("reconstruct fom old state");
#endif

    for (auto i=0; i<n; i++){
      fomStateReconstrObj_(romYprev[i], yFomOld_[i]);
    }

#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("reconstruct fom old state");
#endif
  }


private:
  /* set all entries to zero for all members */
  void resetContainersToZero(){
    yFom_.setZero();
    for (size_t i=0; i<yFomOld_.size(); i++)
      yFomOld_[i].setZero();
  }

protected:
  mutable fom_state_w_t yFom_                             = {};
  mutable std::array<fom_state_w_t, maxNstates> yFomOld_  = {};
  const reconstuctor_type & fomStateReconstrObj_	  = {};

};//end class

}}//end namespace rompp::rom
#endif
