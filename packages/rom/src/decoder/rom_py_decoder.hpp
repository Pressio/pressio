
#ifndef ROM_DECODER_ROM_PY_DECODER_HPP_
#define ROM_DECODER_ROM_PY_DECODER_HPP_

#include "./impl/rom_py_decoder.hpp"

namespace pressio{ namespace rom{

template <
  typename matrix_type,
  typename rom_state_type,
  typename fom_state_type,
  typename ... Args>
using PyDecoder = impl::PyDecoder<matrix_type, rom_state_type, fom_state_type, Args...>;

}} // end namespace pressio::rom
#endif  // ROM_DECODER_ROM_PY_DECODER_HPP_
