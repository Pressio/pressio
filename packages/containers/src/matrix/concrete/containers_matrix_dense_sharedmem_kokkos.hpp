
#ifdef HAVE_KOKKOS
#ifndef CONTAINERS_DENSE_MATRIX_SHAREDMEM_KOKKOS_HPP_
#define CONTAINERS_DENSE_MATRIX_SHAREDMEM_KOKKOS_HPP_

#include "../../shared_base/containers_container_base.hpp"
#include "../base/containers_matrix_sharedmem_base.hpp"

namespace pressio{ namespace containers{

template <typename wrapped_type>
class Matrix<
  wrapped_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_dense_matrix_kokkos<wrapped_type>::value
    >
  >
  : public ContainerBase< Matrix<wrapped_type>, wrapped_type >,
    public MatrixSharedMemBase< Matrix<wrapped_type> >
{

  using this_t = Matrix<wrapped_type>;
  using mytraits = typename details::traits<this_t>;
  using sc_t = typename mytraits::scalar_t;
  using ord_t = typename  mytraits::ordinal_t;
  using wrap_t = typename mytraits::wrapped_t;

  // Views have "view semantics." copy constructor and
  // operator= only do shallow copies.
  // Here, for the time being, we construct wrapper
  // of a view WITHOUT doing shallow copy.
  // We create a new object and deep_copy original.

public:
  Matrix() = default;

  explicit Matrix(const wrap_t src)
    : data_{src.label(),
	    src.extent(0),
	    src.extent(1)}{
    Kokkos::deep_copy(data_, src);
  }

  Matrix(const std::string & label, size_t e1, size_t e2)
    : data_{label, e1, e2}
  {}

  Matrix(const this_t & other)
    : data_{other.data_.label(),
	    other.data_.extent(0),
	    other.data_.extent(1)}{
    Kokkos::deep_copy(data_, other.data_);
  }

  ~Matrix(){}

public:
  // copy assign implments copy semantics not view (for time being)
  this_t & operator=(const this_t & other){
    assert(this->rows() == other.rows());
    assert(this->cols() == other.cols());
    Kokkos::deep_copy(data_, *other.data());
    return *this;
  }

private:
  wrap_t const * dataImpl() const{
    return &data_;
  }
  wrap_t * dataImpl(){
    return &data_;
  }

  wrap_t dataCpImpl(){
    return data_;
  }

  ord_t rowsImpl() const{
    return data_.extent(0);
  }

  ord_t colsImpl() const{
    return data_.extent(1);
  }

private:
  friend ContainerBase< this_t, wrapped_type >;
  friend MatrixSharedMemBase< this_t >;

private:
  wrap_t data_ = {};

};//end class

}}//end namespace pressio::containers
#endif
#endif
