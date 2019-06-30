
#ifndef ROMPP_MPL_IS_SAME_HPP_
#define ROMPP_MPL_IS_SAME_HPP_

namespace rompp { namespace mpl {

template<typename T1, typename T2>
struct is_same : std::is_same<T1,T2> {};

}} // end namespace rompp::mpl

#endif
