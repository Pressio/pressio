
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
	     std::enable_if<
	       core::meta::is_matrixDenseSerialStdlib<wrapped_type
					   >::value
	       >::type
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

public:
  matrix() = delete;
  matrix(ord_t rows, ord_t cols){
    this->resize(rows,cols);
  }
  ~matrix() = default;

public:
  sc_t & operator() (ord_t row, ord_t col){
    // check if we are withinbound 
    return data_(row,col);
  }

  sc_t const & operator() (ord_t row, ord_t col) const{
    // check if we are withinbound 
    return data_(row,col);
  }

private:

  wrap_t const * dataImpl() const{
    return &data_;
  };

  wrap_t * data(){
    return &data_;
  };

  ord_t rows() const{
    assert(!data_.empty());
    return data_.size();
  }

  ord_t cols() const{
    return data_.empty() ? 0 : data_[0].size();
  }

  void resizeImpl(ord_t nrows, ord_t ncols){    
    data_.resize(nrows);
    for (auto & it : data_)
      it.resize(ncols);
  }

private:
  friend matrixGenericBase< matrix<wrapped_type> >;
  friend matrixDenseSerialBase< matrix<wrapped_type> >;
  
private:
  std::vector<std::vector<sc_t>> data_;
 
};
  
}//end namespace core
#endif
