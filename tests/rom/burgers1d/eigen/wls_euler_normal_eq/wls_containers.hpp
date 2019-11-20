namespace pressio{ namespace rom{ namespace wls{
// container for the fom states for WLS
template <std::size_t n, typename fom_state_type, typename reconstuctor_type>
class WlsFomStatesContainer{

  // put here usual things and overload operator [] so we can access the fom state at a given index
  // where the index typically is the step number
public:

  static constexpr std::size_t n_ = n;

  WlsFomStatesContainer() = delete;
  ~WlsFomStatesContainer() = default;

  template <
    typename _fom_state_type = fom_state_type,
    ::pressio::mpl::enable_if_t<
      ::pressio::containers::meta::is_vector_wrapper<_fom_state_type>::value
      > * = nullptr
    >
  WlsFomStatesContainer(const _fom_state_type & fomStateIn,
                        const reconstuctor_type & fomStateReconstr)
    : fomStates_{fomStateIn,fomStateIn}, // does FOM states now hold all fomStateIn?
      fomStateReconstrObj_(fomStateReconstr)
  {
    this->resetContainersToZero();
  }

public:
  fom_state_type & operator[](std::size_t index) {
    return fomStates_[index];
  }

  fom_state_type const & operator[](std::size_t index) const{
    return fomStates_[index];
  }

  template <typename rom_state_t>
  void reconstructFomStateAt(const rom_state_t & romY, std::size_t index) const
  { 
    fomStateReconstrObj_(romY, fomStates_[index]);
  }
  mutable std::array<fom_state_type,n> fomStates_;
  //fom_state_type std::array 
private:
  /* set all entries to zero for all members */
  void resetContainersToZero(){
    for (auto i=0; i<n; i++)
      ::pressio::containers::ops::set_zero(fomStates_[i]);
  }

private:
  const reconstuctor_type & fomStateReconstrObj_  = {};
};
}}}
