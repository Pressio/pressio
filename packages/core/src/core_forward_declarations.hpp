
#ifndef CORE_FORWARDDECLARATIONS_HPP_
#define CORE_FORWARDDECLARATIONS_HPP_

#include "core_ConfigDefs.hpp"

namespace core {


//***************************************
// forward declaration of vector class
//***************************************
template <typename wrapped_type,
	  typename Enable = void>
class vector;

//***************************************
// forward declaration of matrix class
//***************************************
template <typename wrapped_type,
	  typename Enable = void>
class matrix;
  

template <typename vector_a_type,
      typename vector_b_type,
      typename Enable = void>
class dotProduct;


template <typename matrix_type,
      typename vector_type,
      typename Enable = void>
class matVecProduct;

  
} // end namespace core
#endif
