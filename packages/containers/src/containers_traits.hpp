
#ifndef CONTAINERS_TRAITS_HPP_
#define CONTAINERS_TRAITS_HPP_

#include "containers_shared_traits.hpp"

namespace pressio{ namespace containers{ namespace details {

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
