
#ifndef ROMPP_MPL_IS_DEFAULT_CONSTRUCTIBLE_HPP_
#define ROMPP_MPL_IS_DEFAULT_CONSTRUCTIBLE_HPP_

namespace rompp { namespace mpl {

template<typename T>
struct is_default_constructible
: std::is_default_constructible<T> {};

}} // end namespace rompp::mpl

#endif /* ROMPP_MPL_IS_DEFAULT_CONSTRUCTIBLE_HPP_ */
