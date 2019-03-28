
#ifndef ROM_RECONSTRUCTOR_FOM_STATE_HPP_
#define ROM_RECONSTRUCTOR_FOM_STATE_HPP_

#include "rom_ConfigDefs.hpp"

namespace rompp{ namespace rom{

template <typename fom_state_w_t, typename decoder_type>
struct FomStateReconstructor{

  FomStateReconstructor() = delete;

  FomStateReconstructor(const fom_state_w_t & yFomIn,
			const decoder_type & decoder)
    : yFomReference_(yFomIn),
      decoderObj_(decoder){}

  ~FomStateReconstructor() = default;

  template <typename rom_state_t>
  void operator()(const rom_state_t	& romY,
		  fom_state_w_t		& yOut) const
  {
    decoderObj_.applyMapping(romY, yOut);
    yOut += yFomReference_;
  }

  template <typename rom_state_t>
  fom_state_w_t operator()(const rom_state_t & romY) const
  {
    auto yOut(yFomReference_);
    yOut.setZero();
    this->template operator()(romY,yOut);
    return yOut;
  }

private:
  const fom_state_w_t & yFomReference_	= {};
  const decoder_type & decoderObj_	= {};

};//end class

}}//end namespace rompp::rom
#endif
