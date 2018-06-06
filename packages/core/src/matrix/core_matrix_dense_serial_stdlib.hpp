
#ifndef CORE_MATRIX_DENSE_SERIAL_STDLIB_HPP_
#define CORE_MATRIX_DENSE_SERIAL_STDLIB_HPP_

#include <vector>
#include "./base/core_matrix_generic_base.hpp"
#include "./base/core_matrix_dense_serial_base.hpp"



// !!!!!!!!!!!!!!!!!!!!!!!!!!
// WIPL to finish
// !!!!!!!!!!!!!!!!!!!!!!!!!!

namespace core{

template <typename wrapped_type>
class matrix<wrapped_type,
	     typename
	     std::enable_if< core::meta::is_stdlibMatrix<wrapped_type>::value >::type
	     >
  : public matrixGenericBase< matrix<wrapped_type> >,
    public matrixDenseSerialBase< matrix<wrapped_type> >  
{
public:
  using derived_t = matrix<wrapped_type>;
  using sc_t = typename details::traits<derived_t>::scalar_t;
  using ord_t = typename details::traits<derived_t>::ordinal_t;
  using wrap_t = typename details::traits<derived_t>::wrapped_t;
  using der_t = typename details::traits<derived_t>::derived_t;
  
private:
  std::vector<std::vector<sc_t>> data_;

public:
  matrix(){}
  ~matrix(){}
  //-----------------------------------

  wrap_t const * dataImpl() const{
    return &data_;
  };

  wrap_t * data(){
    return &data_;
  };

  size_t rows() const{
    return data_.size();
  }

  size_t cols() const{
    return data_[0].size();
  }

  void resizeImpl(size_t nrows, size_t ncols){    
    data_.resize(nrows);
    for (auto & it : data_)
      it.resize(ncols);
  }
  
};

  
}//end namespace core
#endif
