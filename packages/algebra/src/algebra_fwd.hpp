
#ifndef ALGEBRA_FORWARD_DECLARATIONS_HPP_
#define ALGEBRA_FORWARD_DECLARATIONS_HPP_

#include "algebra_ConfigDefs.hpp"

namespace rompp{ namespace algebra{

template <typename wrapped_type,
	  typename Enable = void>
class Vector;

template <typename wrapped_type,
	  typename Enable = void>
class MultiVector;

template <typename wrapped_type,
	  typename Enable = void>
class Matrix;

template<typename T,
	 typename enable = void>
struct traits;

}} // end namespace rompp::algebra
#endif
