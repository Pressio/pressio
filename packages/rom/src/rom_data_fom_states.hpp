
#ifndef ROM_FOM_STATES_DATA_HPP_
#define ROM_FOM_STATES_DATA_HPP_

#include "rom_ConfigDefs.hpp"

namespace rompp{ namespace rom{

template <typename fom_state_w_type, int maxNstates, typename decoder_type>
struct FomStatesData{

  static constexpr int maxNstates_ = maxNstates;
  using fom_state_w_t = fom_state_w_type;

  FomStatesData() = delete;

  template <int _maxNstates = maxNstates,
            core::meta::enable_if_t<_maxNstates==1> * = nullptr>
  FomStatesData(const fom_state_w_t & yFomIn,
		const decoder_type & decoder)
    : yFomReference_(yFomIn), yFom_(yFomIn), yFomOld_{yFomIn},
      decoderObj_(decoder){}

  template <int _maxNstates = maxNstates,
            core::meta::enable_if_t<_maxNstates==2> * = nullptr>
  FomStatesData(const fom_state_w_t & yFomIn,
		const decoder_type & decoder)
    : yFomReference_(yFomIn), yFom_(yFomIn), yFomOld_{yFomIn, yFomIn},
      decoderObj_(decoder){}

  ~FomStatesData() = default;

protected:
  template <typename rom_state_t>
  void reconstructCurrentFomState(const rom_state_t & romY) const{
    decoderObj_.applyMapping(romY, yFom_);
    yFom_ += yFomReference_;
  }

  template <int n, typename rom_state_t>
  void reconstructFomOldStates(const std::array<rom_state_t, n> & romYprev) const
  {
    for (auto i=0; i<n; i++){
      decoderObj_.applyMapping(romYprev[i], yFomOld_[i]);
      yFomOld_[i] += yFomReference_;
    }
  }

protected:
  const fom_state_w_t & yFomReference_			  = {};
  mutable fom_state_w_t yFom_                             = {};
  mutable std::array<fom_state_w_t, maxNstates> yFomOld_  = {};
  const decoder_type & decoderObj_			  = {};

};//end class

}}//end namespace rompp::rom
#endif
