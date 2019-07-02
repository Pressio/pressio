
#ifndef CONTAINERS_VECTOR_CONCRETE_VECTOR_ARBITRARY_HPP_
#define CONTAINERS_VECTOR_CONCRETE_VECTOR_ARBITRARY_HPP_

#include "../../shared_base/containers_container_base.hpp"
#include <utility>

namespace pressio{ namespace containers{

template <typename wrapped_type>
class Vector<
  wrapped_type,
  mpl::enable_if_t<
    details::traits<Vector<wrapped_type>>::wrapped_vector_identifier
    == details::WrappedVectorIdentifier::Arbitrary
    >
  >
  : public ContainerBase< Vector<wrapped_type>, wrapped_type >
{
  using this_t = Vector<wrapped_type>;

public:
  Vector() = delete;
  ~Vector() = default;

  template <typename ...Args>
  explicit Vector(Args && ... args)
    : data_( std::forward<Args>(args)... ){}

  explicit Vector(const wrapped_type & vecobj)
    : data_(vecobj){}

  Vector(this_t const & other)
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

private:
  wrapped_type data_ = {};

};//end class

}}//end namespace pressio::containers
#endif
