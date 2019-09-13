
#ifndef PRESSIO_MPL_IS_SAME_HPP_
#define PRESSIO_MPL_IS_SAME_HPP_

namespace pressio { namespace mpl {

template<typename T1, typename T2>
struct is_same : std::is_same<T1,T2> {};

}} // end namespace pressio::mpl

#endif
