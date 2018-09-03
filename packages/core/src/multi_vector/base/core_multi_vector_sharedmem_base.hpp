
#ifndef CORE_MULTIVECTOR_BASE_MULTIVECTOR_SHAREDMEM_BASE_HPP_
#define CORE_MULTIVECTOR_BASE_MULTIVECTOR_SHAREDMEM_BASE_HPP_

#include "../core_multi_vector_traits.hpp"

namespace core{
    
template<typename derived_type>
class MultiVectorSharedMemBase
  : private core::details::CrtpBase<MultiVectorSharedMemBase<derived_type>>
{
  static_assert(details::traits<derived_type>::isSharedMem==1,
  "OOPS: distributed concrete vector inheriting from sharedMem base!");
 
private:
  using sc_t = typename details::traits<derived_type>::scalar_t;
  using ord_t = typename details::traits<derived_type>::ordinal_t;

public:
  // ord_t numVectors() const{
  //   return this->underlying().numVectorsImpl();
  // }
  // ord_t length() const {
  //   return this->underlying().lengthImpl();
  // };
  
private:
  using this_t = MultiVectorSharedMemBase<derived_type>;
  friend derived_type;
  friend core::details::CrtpBase<this_t>;

  MultiVectorSharedMemBase() = default;
  ~MultiVectorSharedMemBase() = default;
    
};//end class
  
} // end namespace core
#endif
