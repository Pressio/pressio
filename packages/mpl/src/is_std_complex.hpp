
#ifndef PRESSIO_MPL_IS_STD_COMPLEX_HPP_
#define PRESSIO_MPL_IS_STD_COMPLEX_HPP_

#include <type_traits>
#include <complex>

namespace pressio{ namespace mpl{ 

template <typename T,
	  typename enable = void>
struct is_std_complex : std::false_type{};

template <typename T>
struct is_std_complex<T, typename
		      ::pressio::mpl::enable_if_t<
			   std::is_same<T,
					std::complex<typename T::value_type
						     >
					>::value
			   >
		      > : std::true_type{};

}} // namespace pressio::mpl
#endif
