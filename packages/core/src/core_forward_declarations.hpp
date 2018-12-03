
#ifndef CORE_FORWARD_DECLARATIONS_HPP_
#define CORE_FORWARD_DECLARATIONS_HPP_

#include "core_ConfigDefs.hpp"

namespace rompp{ namespace core{

template <typename wrapped_type,
	  typename Enable = void>
class Vector;

template <typename wrapped_type,
	  typename Enable = void>
class MultiVector;

template <typename wrapped_type,
	  typename Enable = void>
class Matrix;


}} // end namespace rompp::core
#endif
