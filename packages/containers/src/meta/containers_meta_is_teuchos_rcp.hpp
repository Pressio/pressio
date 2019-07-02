
#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_META_IS_TEUCHOS_RCP_HPP_
#define CONTAINERS_META_IS_TEUCHOS_RCP_HPP_

#include <type_traits>
#include <Teuchos_RCPDecl.hpp>

namespace pressio{ namespace containers{ namespace meta {

template <typename T,
	  typename enable = void>
struct is_teuchos_rcp : std::false_type{};

template <typename T>
struct is_teuchos_rcp<
  T, typename
  std::enable_if<
       std::is_same<T,
		    Teuchos::RCP<typename T::element_type>
		    >::value or
       std::is_same<T,
		    Teuchos::RCP<const typename T::element_type>
		    >::value
       >::type
  > : std::true_type{};

}}} // namespace pressio::containers::meta
#endif
#endif
