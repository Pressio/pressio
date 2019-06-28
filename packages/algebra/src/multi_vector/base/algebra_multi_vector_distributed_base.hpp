
#ifndef ALGEBRA_MULTIVECTOR_BASE_MULTIVECTOR_DISTRIBUTED_BASE_HPP_
#define ALGEBRA_MULTIVECTOR_BASE_MULTIVECTOR_DISTRIBUTED_BASE_HPP_

#include "../algebra_multi_vector_traits.hpp"

namespace rompp{ namespace algebra{

template<typename derived_type>
class MultiVectorDistributedBase
  : private utils::details::CrtpBase<
  MultiVectorDistributedBase<derived_type>>
{

  static_assert( details::traits<derived_type>::is_shared_mem==0,
  "OOPS: non-distributed concrete class inheriting \
from multi_vector distributed base!");

private:
  using sc_t = typename details::traits<derived_type>::scalar_t;
  using LO_t = typename details::traits<derived_type>::local_ordinal_t;
  using GO_t = typename details::traits<derived_type>::global_ordinal_t;

public:
  GO_t globalNumVectors() const{
    return this->underlying().globalNumVectorsImpl();
  }

  LO_t localNumVectors() const{
    return this->underlying().localNumVectorsImpl();
  }

  GO_t globalLength() const {
    return this->underlying().globalLengthImpl();
  };

  LO_t localLength() const {
    return this->underlying().localLengthImpl();
  };

  void replaceGlobalValue(GO_t globalRowIndex,
			  GO_t vectorIndex,
			  sc_t value){
    this->underlying().replaceGlobalValueImpl(globalRowIndex, vectorIndex, value);
  }

private:
  /* workaround for nvcc issue with templates, see https://devtalk.nvidia.com/default/topic/1037721/nvcc-compilation-error-with-template-parameter-as-a-friend-within-a-namespace/ */
  template<typename DummyType> struct dummy{using type = DummyType;};
  friend typename dummy<derived_type>::type;

  using this_t = MultiVectorDistributedBase<derived_type>;
  friend utils::details::CrtpBase<this_t>;

  MultiVectorDistributedBase() = default;
  ~MultiVectorDistributedBase() = default;

};//end class

}}//end namespace rompp::algebra
#endif
