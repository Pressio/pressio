
#ifndef PRESSIO_MPL_IS_STD_SHARED_PTR_HPP_
#define PRESSIO_MPL_IS_STD_SHARED_PTR_HPP_

#include <type_traits>
#include <memory>

namespace pressio{ namespace mpl{ 

template <typename T,
	  typename enable = void>
struct is_std_shared_ptr : std::false_type{};

template <typename T>
struct is_std_shared_ptr<
  T, typename
  std::enable_if<
       std::is_same<T,
		    std::shared_ptr<typename T::element_type>
		    >::value or
       std::is_same<T,
		    std::shared_ptr<const typename T::element_type>
		    >::value
       >::type
  > : std::true_type{};

}} // namespace pressio::mpl
#endif
