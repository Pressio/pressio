
#ifndef CONTAINERS_MATRIX_CONCRETE_MATRIX_ARBITRARY_HPP_
#define CONTAINERS_MATRIX_CONCRETE_MATRIX_ARBITRARY_HPP_

#include "../../shared_base/containers_container_base.hpp"

namespace pressio{ namespace containers{

template <typename wrapped_type>
class Matrix<
  wrapped_type,
  mpl::enable_if_t<
    details::traits<Matrix<wrapped_type>>::wrapped_matrix_identifier
    == details::WrappedMatrixIdentifier::Arbitrary
    >
  >
  : public ContainerBase< Matrix<wrapped_type>, wrapped_type >
{

  using this_t = Matrix<wrapped_type>;

public:
  Matrix() = delete;
  ~Matrix() = delete;

  template <typename ...Args>
  Matrix(Args && ... args)
    : data_( std::forward<Args>(args)... ){}

  explicit Matrix(const wrap_t & vecobj)
    : data_(vecobj){}

  Matrix(this_t const & other)
    : data_(*other.data()){}

private:
  wrapped_type const * dataImpl() const{
    return &data_;
  }

  wrapped_type * dataImpl(){
    return &data_;
  }

private:
  friend ContainerBase< this_t, wrapped_type >;

  wrap_t data_ = {};

};//end class

}}//end namespace pressio::containers

#endif
