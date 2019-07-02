
#ifndef CONTAINERS_FORWARD_DECLARATIONS_HPP_
#define CONTAINERS_FORWARD_DECLARATIONS_HPP_

#include "containers_ConfigDefs.hpp"

namespace pressio{ namespace containers{

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

}} // end namespace pressio::containers
#endif
