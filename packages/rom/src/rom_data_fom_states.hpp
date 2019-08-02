
#ifndef ROM_DATA_FOM_STATES_HPP_
#define ROM_DATA_FOM_STATES_HPP_

#include "rom_ConfigDefs.hpp"

namespace pressio{ namespace rom{

template <
  typename fom_state_type,
  int maxNstates,
  typename reconstuctor_type,
  typename enable = void
  >
struct FomStatesData;


template <
  typename fom_state_type,
  int maxNstates,
  typename reconstuctor_type
  >
struct FomStatesData<
  fom_state_type, maxNstates, reconstuctor_type,
  mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper<fom_state_type>::value
    >
  >
{

  static constexpr int maxNstates_ = maxNstates;
  using fom_state_t = fom_state_type;

  FomStatesData() = delete;
  ~FomStatesData() = default;

  template <
    int _maxNstates = maxNstates,
    ::pressio::mpl::enable_if_t<_maxNstates==0> * = nullptr
    >
  FomStatesData(const fom_state_t & yFomIn,
		const reconstuctor_type & fomStateReconstr)
    : yFom_(yFomIn),
      fomStateReconstrObj_(fomStateReconstr)
  {
    this->resetContainersToZero();
  }

  template <
    int _maxNstates = maxNstates,
    ::pressio::mpl::enable_if_t<_maxNstates==1> * = nullptr
    >
  FomStatesData(const fom_state_t & yFomIn,
		const reconstuctor_type & fomStateReconstr)
    : yFom_(yFomIn),
      yFomOld_{{yFomIn}},
      fomStateReconstrObj_(fomStateReconstr)
  {
    this->resetContainersToZero();
  }

  template <
    int _maxNstates = maxNstates,
    ::pressio::mpl::enable_if_t<_maxNstates==2> * = nullptr
    >
  FomStatesData(const fom_state_t & yFomIn,
		const reconstuctor_type & fomStateReconstr)
    : yFom_(yFomIn),
      yFomOld_{{yFomIn, yFomIn}},
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
  mutable fom_state_t yFom_                             = {};
  mutable std::array<fom_state_t, maxNstates> yFomOld_  = {};
  const reconstuctor_type & fomStateReconstrObj_	= {};

};//end class










#ifdef HAVE_PYBIND11
/*
 * partially specialize for pybind array_t
 */
template <
  typename fom_state_type,
  int maxNstates,
  typename reconstuctor_type
  >
struct FomStatesData<
  fom_state_type, maxNstates, reconstuctor_type,
  mpl::enable_if_t<
    ::pressio::containers::meta::is_array_pybind11<fom_state_type>::value
    >
  >
{
  static constexpr int maxNstates_ = maxNstates;
  using fom_state_t = fom_state_type;

  FomStatesData() = delete;
  ~FomStatesData() = default;

  template <
    int _maxNstates = maxNstates,
    ::pressio::mpl::enable_if_t<_maxNstates==0> * = nullptr
    >
  FomStatesData(const fom_state_t & yFomIn,
		const reconstuctor_type & fomStateReconstr)
    : yFom_{ {fom_state_t(const_cast<fom_state_t &>(yFomIn).request())} },
      fomStateReconstrObj_(fomStateReconstr)
  {
    this->resetContainersToZero();
  }

  template <
    int _maxNstates = maxNstates,
    ::pressio::mpl::enable_if_t<_maxNstates==1> * = nullptr
    >
  FomStatesData(const fom_state_t & yFomIn,
		const reconstuctor_type & fomStateReconstr)
    : yFom_{ {fom_state_t(const_cast<fom_state_t &>(yFomIn).request())} },
      yFomOld_{ {fom_state_t(const_cast<fom_state_t &>(yFomIn).request())} },
      fomStateReconstrObj_(fomStateReconstr)
  {
    this->resetContainersToZero();

#ifdef DEBUG_PRINT
    std::cout << "inside fom states" << std::endl;
    std::cout << "FomStatesData: yFom_ " << yFom_.data() << std::endl;
    std::cout << "FomStatesData: yFomOld[0]_ " << yFomOld_[0].data() << std::endl;
    std::cout << std::endl;
#endif
  }

  template <
    int _maxNstates = maxNstates,
    ::pressio::mpl::enable_if_t<_maxNstates==2> * = nullptr
    >
  FomStatesData(const fom_state_t & yFomIn,
		const reconstuctor_type & fomStateReconstr)
    : yFom_{ {fom_state_t(const_cast<fom_state_t &>(yFomIn).request())} },
      yFomOld_{ {fom_state_t(const_cast<fom_state_t &>(yFomIn).request())},
		{fom_state_t(const_cast<fom_state_t &>(yFomIn).request())}},
      fomStateReconstrObj_(fomStateReconstr)
  {
    this->resetContainersToZero();
  }

protected:

  template <typename rom_state_t>
  void reconstructCurrentFomState(const rom_state_t & romY) const{
    fomStateReconstrObj_(romY, yFom_);
  }

  template <int n, typename rom_state_t>
  void reconstructFomOldStates(const std::array<
			       rom_state_t, n
			       > & romYprev) const{
    for (auto i=0; i<n; i++){
      fomStateReconstrObj_(romYprev[i], yFomOld_[i]);
    }
  }

private:
  /* set all entries to zero for all members */
  void resetContainersToZero(){
    ::pressio::containers::ops::set_zero(yFom_);
    for (size_t i=0; i<yFomOld_.size(); i++)
      ::pressio::containers::ops::set_zero(yFomOld_[i]);
  }

protected:
  mutable fom_state_t yFom_                             = {};
  mutable std::array<fom_state_t, maxNstates> yFomOld_  = {};
  const reconstuctor_type & fomStateReconstrObj_	= {};

};//end class
#endif


}}//end namespace pressio::rom
#endif
