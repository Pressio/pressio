
#ifndef ALGEBRA_META_META_BASIC_HPP_
#define ALGEBRA_META_META_BASIC_HPP_

#include "algebra_meta_has_communicator_typedef.hpp"
#include "algebra_meta_has_data_map_typedef.hpp"
#include "algebra_meta_has_global_ordinal_typedef.hpp"
#include "algebra_meta_has_local_ordinal_typedef.hpp"
#include "algebra_meta_has_ordinal_typedef.hpp"
#include "algebra_meta_has_scalar_typedef.hpp"
#include "algebra_meta_has_size_method.hpp"
#include "algebra_meta_is_teuchos_rcp.hpp"


namespace rompp{ namespace algebra{ namespace meta {

template<typename T>
struct remove_const: std::remove_const<T>{};

template<typename T>
struct remove_reference : std::remove_reference<T>{};

template<typename T>
struct remove_pointer : std::remove_pointer<T>{};

template<typename T>
struct is_arithmetic : std::is_arithmetic<T>{};

template<typename T>
struct is_integral: std::is_integral<T>{};

}}} // namespace rompp::algebra::meta
#endif
