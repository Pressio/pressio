
#ifndef CORE_VECTOR_BASE_VECTOR_SHAREDMEM_BASE_HPP_
#define CORE_VECTOR_BASE_VECTOR_SHAREDMEM_BASE_HPP_

#include "../core_vector_traits.hpp"

namespace core{
    
template<typename derived_type>
class VectorSharedMemBase
  : private core::details::CrtpBase<VectorSharedMemBase<derived_type>>
{

  static_assert(details::traits<derived_type>::isSharedMem==1,
  "OOPS: distributed concrete vector inheriting from sharedMem base!");
  
private:
  using this_t = VectorSharedMemBase<derived_type>;
  using ord_t = typename details::traits<derived_type>::ordinal_t;

public:
  ord_t size() const {
    return this->underlying().sizeImpl();
  };
  void resize(ord_t newSize) {
    this->underlying().resizeImpl(newSize);
  };
    
private:
  friend derived_type;
  friend core::details::CrtpBase<this_t>;

  VectorSharedMemBase() = default;
  ~VectorSharedMemBase() = default;
    
};//end class
  
} // end namespace core
#endif
