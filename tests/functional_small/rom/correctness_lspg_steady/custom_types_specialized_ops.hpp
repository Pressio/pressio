
#ifndef PRESSIO_TESTS_LSPG_STEADY_CORRECT_SPECIALIZED_OPS_CUSTOM_TYPES_HPP_
#define PRESSIO_TESTS_LSPG_STEADY_CORRECT_SPECIALIZED_OPS_CUSTOM_TYPES_HPP_

#include "../custom_data_types.hpp"
#include "pressio/ops.hpp"

namespace pressio{ 

template<class ScalarType> 
struct Traits<pressiotests::MyCustomVector<ScalarType>>{
  using scalar_type = ScalarType;
};

template<class ScalarType> 
struct Traits<pressiotests::MyCustomMatrix<ScalarType>>{
  using scalar_type = ScalarType;
};

namespace ops{

template<class ScalarType> 
std::size_t extent(pressiotests::MyCustomVector<ScalarType> & object, int i){
    return object.extent(i);
}

template<class ScalarType> 
std::size_t extent(pressiotests::MyCustomMatrix<ScalarType> & object, int i){
    return object.extent(i);
}

template<class ScalarType> 
void set_zero(pressiotests::MyCustomVector<ScalarType> & object){
  object.fill(0);
}

template<class ScalarType> 
void set_zero(pressiotests::MyCustomMatrix<ScalarType> & object){
  object.fill(0);
}

template<class ScalarType> 
void deep_copy(pressiotests::MyCustomVector<ScalarType> & dest, 
               const pressiotests::MyCustomVector<ScalarType> & src){
  dest = src;
}

template<class ScalarType> 
void deep_copy(pressiotests::MyCustomMatrix<ScalarType> & dest, 
               const pressiotests::MyCustomMatrix<ScalarType> & src){
  dest = src;
}

template<class ScalarType> 
pressiotests::MyCustomVector<ScalarType> clone(const pressiotests::MyCustomVector<ScalarType> & src){
  return pressiotests::MyCustomVector<ScalarType>(src.extent(0));
}

template<class ScalarType> 
pressiotests::MyCustomMatrix<ScalarType> clone(const pressiotests::MyCustomMatrix<ScalarType> & src){
  return pressiotests::MyCustomMatrix<ScalarType>(src.extent(0), src.extent(1));
}

// this update is needed by the decoder
template<class ScalarType> 
void update(pressiotests::MyCustomVector<ScalarType> & v,        const ScalarType a,
            const pressiotests::MyCustomVector<ScalarType> & v1, const ScalarType b)
{
  for (std::size_t i=0; i< v.extent(0); ++i){
    v(i) = a*v(i) + b*v1(i);
  }
}

// this product is needed by the decoder
// z = beta*z + alpha * A * x
// where x is indexable as x(i)
template<class x_t, class ScalarType>
void product(pressio::nontranspose,
             ScalarType alpha,
             const pressiotests::MyCustomMatrix<ScalarType> & A,
             const x_t & x,
             ScalarType beta,
             pressiotests::MyCustomVector<ScalarType> & z)
{
  // obviously not efficient, just for demonstration
  for (std::size_t i=0; i<A.extent(0); ++i)
  {
    z(i) = beta*z(i);
    for (std::size_t j=0; j<A.extent(1); ++j){
      z(i) += alpha*A(i,j)*x(j);
    }
  }
}

}}//end namespace pressio::ops

#endif
