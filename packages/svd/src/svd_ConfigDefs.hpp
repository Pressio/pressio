
#ifndef SVD_CONFIGDEFS_HPP
#define SVD_CONFIGDEFS_HPP

#include "core_ConfigDefs.hpp"
#include "svd_config.h"

namespace ode{
namespace details {

template<typename T, typename enable = void>
struct traits : core::details::traits<T> {};
  
} // end namespace details
  
} // end of svd namespace

#endif
