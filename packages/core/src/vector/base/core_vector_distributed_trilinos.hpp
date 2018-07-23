
#ifndef CORE_VECTOR_BASE_VECTOR_DISTRIBUTED_TRILINOS_HPP_
#define CORE_VECTOR_BASE_VECTOR_DISTRIBUTED_TRILINOS_HPP_

#include "../core_vector_traits.hpp"

namespace core{
    
template<typename derived_type>
class VectorDistributedTrilinos
  : private core::details::CrtpBase<VectorDistributedTrilinos<derived_type>>
{

  static_assert( details::traits<derived_type>::isDistributed==1,
  "OOPS: serial concrete vector inheriting from distributed trilonos!");

private:
  using this_t = VectorDistributedTrilinos<derived_type>;
  using map_t = typename details::traits<derived_type>::data_map_t;
    
public:
  map_t const & getDataMap() const{
    return this->underlying().getDataMapImpl();
  }

  void replaceDataMap(const map_t & mapObj){
    return this->underlying().replaceDataMapImpl(mapObj);
  }
  
private:
  friend derived_type;
  friend core::details::CrtpBase<this_t>;

  VectorDistributedTrilinos() = default;
  ~VectorDistributedTrilinos() = default;

};//end class
} // end namespace core
#endif
