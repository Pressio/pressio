
#ifdef DEBUG_PRINT
#ifndef CORE_SET_STREAM_PRECISION_HPP_
#define CORE_SET_STREAM_PRECISION_HPP_

#include <iomanip>

namespace rompp{ namespace core{ namespace impl{

template <typename stream_t, typename scalar_t>
void setStreamPrecision(stream_t & ss){
  constexpr auto prec = std::is_same<scalar_t,
				     double>::value ? 15 : 6;
  ss << std::setprecision(prec);
}

}}} // end of namespace rompp::core::impl
#endif
#endif

