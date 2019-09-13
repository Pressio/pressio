
#ifndef PRESSIO_MPL_NOT_SAME_HPP_
#define PRESSIO_MPL_NOT_SAME_HPP_

namespace pressio { namespace mpl {

template<class T, class U>
struct not_same : std::true_type {};

template<class T>
struct not_same<T, T> : std::false_type {};

}} // end namespace pressio::mpl

#endif
