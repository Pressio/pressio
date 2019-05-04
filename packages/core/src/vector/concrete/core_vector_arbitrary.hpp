
#ifndef CORE_VECTOR_CONCRETE_VECTOR_ARBITRARY_HPP_
#define CORE_VECTOR_CONCRETE_VECTOR_ARBITRARY_HPP_

#include "../../shared_base/core_container_base.hpp"
#include <utility>

namespace rompp{ namespace core{

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
  using wrap_t = typename details::traits<this_t>::wrapped_t;

public:
  Vector() = delete;
  ~Vector() = default;

  template <typename ...Args>
  Vector(Args && ... args)
    : data_( std::forward<Args>(args)... ){}

  explicit Vector(const wrap_t & vecobj)
    : data_(vecobj){}

  Vector(this_t const & other)
    : data_(*other.data()){}

private:

  wrap_t const * dataImpl() const{
    return &data_;
  }

  wrap_t * dataImpl(){
    return &data_;
  }

private:
  friend ContainerBase< this_t, wrapped_type >;

private:
  wrap_t data_ = {};

};//end class

}}//end namespace rompp::core
#endif
