
#ifndef CORE_MATRIX_CONCRETE_MATRIX_DENSE_SHAREDMEM_STDLIB_HPP_
#define CORE_MATRIX_CONCRETE_MATRIX_DENSE_SHAREDMEM_STDLIB_HPP_

#include "../../shared_base/core_container_base.hpp"
#include "../base/core_matrix_sharedmem_base.hpp"
#include "../base/core_matrix_dense_sharedmem_base.hpp"


// !!!!!!!!!!!!!!!!!!!!!!!!!!
// WIP to finish
// !!!!!!!!!!!!!!!!!!!!!!!!!!

namespace core{

template <typename wrapped_type>
class Matrix<wrapped_type,
       typename
       std::enable_if<
	 core::meta::is_matrix_dense_sharedmem_stdlib<
	   wrapped_type>::value
	 >::type
       >
  : public ContainerBase< Matrix<wrapped_type>, wrapped_type >,
    public MatrixDenseSharedMemBase< Matrix<wrapped_type> >
{

  using derived_t = Matrix<wrapped_type>;
  using sc_t = typename details::traits<derived_t>::scalar_t;
  using ord_t = typename details::traits<derived_t>::ordinal_t;
  using wrap_t = typename details::traits<derived_t>::wrapped_t;
  using der_t = typename details::traits<derived_t>::derived_t;

public:
  Matrix() = delete;
  Matrix(ord_t rows, ord_t cols){
    this->resize(rows,cols);
  }
  ~Matrix() = default;

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

  // !!!!!!!!!!!!!!!!!!!!!!!!!!
  // WIP to finish
  // !!!!!!!!!!!!!!!!!!!!!!!!!!

  
private:
  friend ContainerBase< Matrix<wrapped_type>, wrapped_type >;
  friend MatrixDenseSharedMemBase< Matrix<wrapped_type> >;
  
private:
  std::vector<std::vector<sc_t>> data_;
 
};
  
}//end namespace core
#endif
