
#ifndef CORE_TRAITS_HPP_
#define CORE_TRAITS_HPP_

#include "core_shared_traits.hpp"

namespace rompp{ namespace core{ namespace details {

template<typename T, typename enable = void>
struct traits : public
containers_shared_traits<void, void,
			 false, false, false,
			 WrappedPackageIdentifier::Undefined,
			 false, false>{};

template<typename T>
struct traits<const T> : traits<T> {};

}}} // end namespace
#endif
