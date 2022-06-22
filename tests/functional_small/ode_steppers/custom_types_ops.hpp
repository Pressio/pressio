
#ifndef ODE_TESTING_FUNCTIONAL_SMALL_STEPPERS_CUSTOM_OPS_HPP_
#define ODE_TESTING_FUNCTIONAL_SMALL_STEPPERS_CUSTOM_OPS_HPP_

#include "pressio/type_traits.hpp"

namespace pressio{

template<>
struct Traits<VectorType>
{
  using scalar_type = ScalarType;
};

namespace ops{
void deep_copy(VectorType & dest, const VectorType & from){
  dest = from;
}

VectorType clone(const VectorType & src)
{
  return VectorType(src.size());
}

void update(VectorType & v,        const ScalarType a,
	    const VectorType & v1, const ScalarType b)
{
  for (size_t i=0; i<v.size(); ++i){
    v[i] = a*v[i] + b*v1[i];
  }
}

void update(VectorType & v,        const ScalarType a,
            const VectorType & v1, const ScalarType b,
            const VectorType & v2, const ScalarType c)
{
  for (size_t i=0; i<v.size(); ++i){
    v[i] = a*v[i] + b*v1[i] + c*v2[i];
  }
}

void update(VectorType & v, const ScalarType a,
    const VectorType & v1, const ScalarType b,
    const VectorType & v2, const ScalarType c,
    const VectorType & v3, const ScalarType d,
    const VectorType & v4, const ScalarType e)
{
  for (size_t i=0; i<v.size(); ++i){
    v[i] = a*v[i] + b*v1[i] + c*v2[i] + d*v3[i] + e*v4[i];
  }
}
}} //end namespace pressio::ops

#endif
