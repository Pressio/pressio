
#ifndef CORE_MATRIX_CONCRETE_MATRIX_DENSE_SHAREDMEM_EIGEN_HPP_
#define CORE_MATRIX_CONCRETE_MATRIX_DENSE_SHAREDMEM_EIGEN_HPP_

#include "../../shared_base/core_container_base.hpp"
#include "../../shared_base/core_container_resizable_base.hpp"
#include "../../shared_base/core_container_nonresizable_base.hpp"
#include "../../shared_base/core_container_subscriptable_base.hpp"

#include "../base/core_matrix_base.hpp"
#include "../base/core_matrix_sharedmem_base.hpp"
#include "../base/core_matrix_dense_sharedmem_base.hpp"
#include "../base/core_matrix_math_base.hpp"

namespace rompp{
namespace core{

template <typename wrapped_type>
class Matrix<wrapped_type,
	     core::meta::enable_if_t<
	       core::meta::is_matrix_dense_sharedmem_eigen<
		 wrapped_type>::value>
	     >
  : public ContainerBase< Matrix<wrapped_type>, wrapped_type >,
    public MatrixBase< Matrix<wrapped_type> >,
    public MatrixSharedMemBase< Matrix<wrapped_type> >,
    public MatrixDenseSharedMemBase< Matrix<wrapped_type> >,
    public ContainerSubscriptable2DBase<
     Matrix<wrapped_type>, 
     typename details::traits<Matrix<wrapped_type>>::scalar_t,
     typename details::traits<Matrix<wrapped_type>>::ordinal_t>,
    public std::conditional<
      details::traits<Matrix<wrapped_type>>::is_static == true,
      ContainerNonResizableBase<Matrix<wrapped_type>, 2>,
      ContainerResizableBase<Matrix<wrapped_type>, 2>>::type
{

  using derived_t = Matrix<wrapped_type>;
  using mytraits = typename details::traits<derived_t>;  
  using sc_t = typename mytraits::scalar_t;
  using ord_t = typename mytraits::ordinal_t;
  using wrap_t = typename mytraits::wrapped_t;
  using der_t = typename mytraits::derived_t;
  
public:
  Matrix() = default;

  template <typename T = ord_t,
	    typename std::enable_if<
	      !mytraits::is_static, T
	      >::type * = nullptr>
  explicit Matrix(T nrows, T ncols) {
    this->resize(nrows,ncols);
  }

  explicit Matrix(const wrap_t & other) : data_(other){}

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
  
  derived_t & operator+=(const derived_t & other) {
    assert( other.rows() == this->rows() );
    assert( other.cols() == this->cols() );
    this->data_ += *other.data();
    return *this;
  }

  derived_t & operator-=(const derived_t & other) {
    assert( other.rows() == this->rows() );
    assert( other.cols() == this->cols() );
    this->data_ -= *other.data();
    return *this;
  }

private:  
  wrap_t * dataImpl(){
    return &data_;
  };  

  wrap_t const * dataImpl() const{
    return &data_;
  };

  void addToDiagonalImpl(sc_t value) {
    // check matrix is diagonal
    assert(this->rows()==this->cols());
    for (ord_t ir=0; ir<this->rows(); ir++)
      data_(ir,ir) += value;
  };
  
  ord_t rowsImpl() const{
    return data_.rows();
  }

  ord_t colsImpl() const{
    return data_.cols();
  }

  template <typename T = ord_t,
	    typename std::enable_if<
	      !mytraits::is_static, T
	      >::type * = nullptr>
  void resizeImpl(T nrows, T ncols){
    data_.resize(nrows, ncols);
    data_ = wrapped_type::Zero(nrows,ncols);
  }
    
private:
  friend ContainerBase< derived_t, wrapped_type >;
  friend MatrixBase< derived_t >;
  friend MatrixSharedMemBase< derived_t >;
  friend MatrixDenseSharedMemBase< derived_t >;
  friend ContainerSubscriptable2DBase< derived_t, sc_t, ord_t>;
  friend typename std::conditional<
    details::traits<derived_t>::is_static == true,
    ContainerNonResizableBase<derived_t, 2>,
    ContainerResizableBase<derived_t, 2>
    >::type;

private:
  wrap_t data_;
     
};//end class 
}//end namespace core 
}//end namespace rompp
#endif
