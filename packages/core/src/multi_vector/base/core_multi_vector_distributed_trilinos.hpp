
#ifndef CORE_MULTIVECTOR_BASE_MULTIVECTOR_DISTRIBUTED_TRILINOS_HPP_
#define CORE_MULTIVECTOR_BASE_MULTIVECTOR_DISTRIBUTED_TRILINOS_HPP_

#include "../core_multi_vector_traits.hpp"

namespace core{
    
template<typename derived_type>
class MultiVectorDistributedTrilinos
  : private core::details::CrtpBase<MultiVectorDistributedTrilinos<derived_type>>
{

  static_assert( details::traits<derived_type>::isDistributed==1,
  "OOPS: serial concrete vector inheriting from distributed trilonos!");

private:
  using this_t = MultiVectorDistributedTrilinos<derived_type>;
  using map_t = typename details::traits<derived_type>::data_map_t;
    
// public:
//   map_t const & getDataMap() const{
//     return this->underlying().getDataMapImpl();
//   }

//   void replaceDataMap(const map_t & mapObj){
//     return this->underlying().replaceDataMapImpl(mapObj);
//   }
  
private:
  friend derived_type;
  friend core::details::CrtpBase<this_t>;

  MultiVectorDistributedTrilinos() = default;
  ~MultiVectorDistributedTrilinos() = default;

};//end class
} // end namespace core
#endif
