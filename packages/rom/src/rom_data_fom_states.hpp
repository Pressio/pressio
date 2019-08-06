
#ifndef ROM_DATA_FOM_STATES_HPP_
#define ROM_DATA_FOM_STATES_HPP_

#include "rom_ConfigDefs.hpp"
#include "rom_fwd.hpp"

namespace pressio{ namespace rom{

template <
  typename fom_state_type,
  int maxNstates,
  typename reconstuctor_type
  >
struct FomStatesData<fom_state_type, maxNstates, reconstuctor_type>
{
  static constexpr int maxNstates_ = maxNstates;
  //using fom_state_t = fom_state_type;

  FomStatesData() = delete;
  ~FomStatesData() = default;

  /* ----------------
   * maxNstates = 0
   * ---------------*/

  /* cnstr for maxNstates = 0, and fom_state_type is a pressio vector wrapper */
  template <
    typename _fom_state_type = fom_state_type,
    int _maxNstates = maxNstates,
    ::pressio::mpl::enable_if_t<
      _maxNstates==0 and
      ::pressio::containers::meta::is_vector_wrapper<_fom_state_type>::value
      > * = nullptr
    >
  FomStatesData(const _fom_state_type & yFomIn,
		const reconstuctor_type & fomStateReconstr)
    : yFom_(yFomIn),
      fomStateReconstrObj_(fomStateReconstr)
  {
    this->resetContainersToZero();
  }

#ifdef HAVE_PYBIND11
  /* cnstr for maxNstates = 0, and fom_state_type is a pybind11 array */
  template <
    typename _fom_state_type = fom_state_type,
    int _maxNstates = maxNstates,
    ::pressio::mpl::enable_if_t<
      _maxNstates==0 and
      ::pressio::containers::meta::is_array_pybind11<_fom_state_type>::value
      > * = nullptr
    >
  FomStatesData(const _fom_state_type & yFomIn,
		const reconstuctor_type & fomStateReconstr)
    : yFom_{ {_fom_state_type(const_cast<_fom_state_type &>(yFomIn).request())} },
      fomStateReconstrObj_(fomStateReconstr)
  {
    this->resetContainersToZero();
  }
#endif


  /* ----------------
   * maxNstates = 1
   * ---------------*/

  /* cnstr for maxNstates = 1, and fom_state_type is a pressio vector wrapper */
  template <
    typename _fom_state_type = fom_state_type,
    int _maxNstates = maxNstates,
    ::pressio::mpl::enable_if_t<
      _maxNstates==1 and
      ::pressio::containers::meta::is_vector_wrapper<_fom_state_type>::value
      > * = nullptr
    >
  FomStatesData(const _fom_state_type & yFomIn,
		const reconstuctor_type & fomStateReconstr)
    : yFom_(yFomIn),
      yFomOld_{{yFomIn}},
      fomStateReconstrObj_(fomStateReconstr)
  {
    this->resetContainersToZero();
  }

#ifdef HAVE_PYBIND11
  /* cnstr for maxNstates = 1, and fom_state_type is a pybind11 array */
  template <
    typename _fom_state_type = fom_state_type,
    int _maxNstates = maxNstates,
    ::pressio::mpl::enable_if_t<
      _maxNstates==1 and
      ::pressio::containers::meta::is_array_pybind11<_fom_state_type>::value
      > * = nullptr
    >
  FomStatesData(const _fom_state_type & yFomIn,
		const reconstuctor_type & fomStateReconstr)
    : yFom_{ {_fom_state_type(const_cast<_fom_state_type &>(yFomIn).request())} },
      yFomOld_{ {_fom_state_type(const_cast<_fom_state_type &>(yFomIn).request())} },
      fomStateReconstrObj_(fomStateReconstr)
  {
    this->resetContainersToZero();
  }
#endif



  /* ----------------
   * maxNstates = 2
   * ---------------*/

  /* cnstr for maxNstates = 2, and fom_state_type is a pressio vector wrapper */
  template <
    typename _fom_state_type = fom_state_type,
    int _maxNstates = maxNstates,
    ::pressio::mpl::enable_if_t<
      _maxNstates==2 and
      ::pressio::containers::meta::is_vector_wrapper<_fom_state_type>::value
      > * = nullptr
    >
  FomStatesData(const _fom_state_type & yFomIn,
		const reconstuctor_type & fomStateReconstr)
    : yFom_(yFomIn),
      yFomOld_{{yFomIn, yFomIn}},
      fomStateReconstrObj_(fomStateReconstr)
  {
    this->resetContainersToZero();
  }

#ifdef HAVE_PYBIND11
  /* cnstr for maxNstates = 2, and fom_state_type is a pybind11 array */
  template <
    typename _fom_state_type = fom_state_type,
    int _maxNstates = maxNstates,
    ::pressio::mpl::enable_if_t<
      _maxNstates==2 and
      ::pressio::containers::meta::is_array_pybind11<_fom_state_type>::value
      > * = nullptr
    >
  FomStatesData(const _fom_state_type & yFomIn,
		const reconstuctor_type & fomStateReconstr)
    : yFom_{ {_fom_state_type(const_cast<_fom_state_type &>(yFomIn).request())} },
      yFomOld_{ {_fom_state_type(const_cast<_fom_state_type &>(yFomIn).request())},
		{_fom_state_type(const_cast<_fom_state_type &>(yFomIn).request())}},
      fomStateReconstrObj_(fomStateReconstr)
  {
    this->resetContainersToZero();
  }
#endif


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
    ::pressio::containers::ops::set_zero(yFom_);
    for (size_t i=0; i<yFomOld_.size(); i++)
      ::pressio::containers::ops::set_zero(yFomOld_[i]);
  }

protected:
  mutable fom_state_type yFom_                             = {};
  mutable std::array<fom_state_type, maxNstates> yFomOld_  = {};
  const reconstuctor_type & fomStateReconstrObj_	   = {};

};//end class

}}//end namespace pressio::rom
#endif
