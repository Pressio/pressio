
#ifndef ROM_RECONSTRUCTOR_FOM_STATE_HPP_
#define ROM_RECONSTRUCTOR_FOM_STATE_HPP_

#include "rom_ConfigDefs.hpp"

namespace rompp{ namespace rom{

template <
  typename fom_state_t,
  typename decoder_type,
  typename enable = void
  >
struct FomStateReconstructor;


template <
  typename fom_state_t,
  typename decoder_type
  >
struct FomStateReconstructor<
  fom_state_t, decoder_type,
  mpl::enable_if_t<
    ::rompp::core::meta::is_core_wrapper<fom_state_t>::value
    >
  >
{

  FomStateReconstructor() = delete;

  FomStateReconstructor(const fom_state_t & yFomIn,
			const decoder_type & decoder)
    : yFomReference_(yFomIn),
      decoderObj_(decoder)
  {}

  ~FomStateReconstructor() = default;

  template <typename rom_state_t>
  void operator()(const rom_state_t	& romY,
		  fom_state_t		& yOut) const {
    decoderObj_.applyMapping(romY, yOut);
    yOut += yFomReference_;
  }

  template <typename rom_state_t>
  fom_state_t operator()(const rom_state_t & romY) const{
    auto yOut(yFomReference_);
    yOut.setZero();
    this->template operator()(romY,yOut);
    return yOut;
  }

private:
  const fom_state_t & yFomReference_	= {};
  const decoder_type & decoderObj_	= {};

};//end class






#ifdef HAVE_PYBIND11
template <
  typename fom_state_t,
  typename decoder_type
  >
class FomStateReconstructor<
  fom_state_t, decoder_type,
  mpl::enable_if_t<
    ::rompp::core::meta::is_cstyle_array_pybind11<fom_state_t>::value
    >
  >
{
  using scalar_t = typename fom_state_t::value_type;

public:
  FomStateReconstructor() = delete;
  ~FomStateReconstructor() = default;

  FomStateReconstructor(const fom_state_t & yFomIn,
			const decoder_type & decoder)
    : /*yFomReference_( const_cast<fom_state_t&>(yFomIn).request().shape,
		      const_cast<fom_state_t&>(yFomIn).request().strides,
		      yFomIn.data() ),*/
      yFomReference_(yFomIn),
      decoderObj_(decoder)
  {
    std::cout << std::endl;
    std::cout << "Inside FomStateReconstructor " << std::endl;
    std::cout << "yFomReference_ " << yFomReference_.data() << std::endl;
    std::cout << std::endl;
  }

public:
  template <typename rom_state_t>
  void operator()(const rom_state_t	& romY,
		  fom_state_t		& yOut) const {
    decoderObj_.applyMapping(romY, yOut);

    constexpr auto one = ::rompp::core::constants::one<scalar_t>();
    ::rompp::core::ops::do_update(yOut, one, yFomReference_, one);
    //yOut += yFomReference_;
  }

  template <typename rom_state_t>
  fom_state_t operator()(const rom_state_t & romY) const{
    fom_state_t yOut{ fom_state_t(yFomReference_.request()) };
    ::rompp::core::ops::set_zero(yOut);
    this->template operator()(romY,yOut);
    return yOut;
  }

private:
  fom_state_t yFomReference_	= {};
  const decoder_type & decoderObj_	= {};

};//end class
#endif

}}//end namespace rompp::rom
#endif
