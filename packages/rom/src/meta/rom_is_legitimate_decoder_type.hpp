
#ifndef ROM_IS_LEGITIMATE_DECODER_TYPE_HPP_
#define ROM_IS_LEGITIMATE_DECODER_TYPE_HPP_

#include "../decoder/rom_decoder_base.hpp"

namespace pressio{ namespace rom{ namespace meta {

template<typename T, typename enable = void>
struct is_legitimate_decoder_type
  : std::false_type{};

template <typename T>
struct is_legitimate_decoder_type<
  T,
  ::pressio::mpl::enable_if_t<
    ::pressio::mpl::publicly_inherits_from<
      T, ::pressio::rom::DecoderBase<T, typename T::jacobian_t>
      >::value
    >
  > : std::true_type{};

}}} // namespace pressio::rom::meta
#endif
