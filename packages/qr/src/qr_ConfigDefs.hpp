
#ifndef QR_CONFIGDEFS_HPP
#define QR_CONFIGDEFS_HPP

#include "qr_config.h"
#include "../../containers/src/containers_ConfigDefs.hpp"

namespace rompp{ namespace qr{ namespace details {

template<typename T, typename enable = void> struct traits{};
template<typename T>  struct traits<const T> : traits<T> {};

}}}// end namespace rompp::qr::details
#endif
