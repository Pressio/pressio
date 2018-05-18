
#ifndef CORE_FORWARDDECLARATIONS_HPP
#define CORE_FORWARDDECLARATIONS_HPP


#include "core_ConfigDefs.hpp"


namespace core {
namespace details {

  template<typename T, typename enable = void>
  struct traits;

  // here we say once and for all that traits<const T> == traits<T>
  // When constness must affect traits, has to be constness on templates on which T depends.
  // For example, traits<Map<const T> > != traits<Map<T> >, but
  //              traits<const Map<T> > == traits<Map<T> >
  template<typename T> 
  struct traits<const T> : traits<T> {};

} // end namespace details


  //***************************************
  // forward declaration of vector class
  //***************************************
  template <typename wrapped_type, typename enable = void>
  class vector;
  
  //***************************************
  // forward declaration of matrix class
  //***************************************
  template <typename wrapped_type, typename enable = void>
  class matrix;
  
  
} // end namespace core

#endif
