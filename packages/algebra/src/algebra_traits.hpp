
#ifndef ALGEBRA_TRAITS_HPP_
#define ALGEBRA_TRAITS_HPP_

#include "algebra_shared_traits.hpp"

namespace rompp{ namespace algebra{ namespace details {

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
