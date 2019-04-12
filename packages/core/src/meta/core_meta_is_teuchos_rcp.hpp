
#ifdef HAVE_TRILINOS
#ifndef CORE_META_IS_TEUCHOS_RCP_HPP_
#define CORE_META_IS_TEUCHOS_RCP_HPP_

#include <type_traits>
#include <Teuchos_RCPDecl.hpp>

namespace rompp{ namespace core{ namespace meta {

template <typename T,
	  typename enable = void>
struct is_teuchos_rcp_ptr : std::false_type{};

template <typename T>
struct is_teuchos_rcp_ptr<
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

}}} // namespace rompp::core::meta
#endif
#endif
