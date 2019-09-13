
#ifndef PRESSIO_MPL_IS_DEFAULT_CONSTRUCTIBLE_HPP_
#define PRESSIO_MPL_IS_DEFAULT_CONSTRUCTIBLE_HPP_

namespace pressio { namespace mpl {

template<typename T>
struct is_default_constructible
: std::is_default_constructible<T> {};

}} // end namespace pressio::mpl

#endif /* PRESSIO_MPL_IS_DEFAULT_CONSTRUCTIBLE_HPP_ */
