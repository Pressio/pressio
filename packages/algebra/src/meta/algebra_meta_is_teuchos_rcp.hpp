
#ifdef HAVE_TRILINOS
#ifndef ALGEBRA_META_IS_TEUCHOS_RCP_HPP_
#define ALGEBRA_META_IS_TEUCHOS_RCP_HPP_

#include <type_traits>
#include <Teuchos_RCPDecl.hpp>

namespace rompp{ namespace algebra{ namespace meta {

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

}}} // namespace rompp::algebra::meta
#endif
#endif
