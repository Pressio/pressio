
#ifndef CORE_MULTIVECTOR_CONCRETE_MULTIVECTOR_ARBITRARYHPP_
#define CORE_MULTIVECTOR_CONCRETE_MULTIVECTOR_ARBITRARYHPP_

#include "../../shared_base/core_container_base.hpp"

namespace rompp{ namespace core{

template <typename wrapped_type>
class MultiVector<
  wrapped_type,
  mpl::enable_if_t<
    details::traits<MultiVector<wrapped_type>>::wrapped_multi_vector_identifier
    == details::WrappedMultiVectorIdentifier::Arbitrary
    >
  >
  : public ContainerBase< MultiVector<wrapped_type>, wrapped_type >
{

  using this_t = MultiVector<wrapped_type>;

public:
  MultiVector() = delete;
  ~MultiVector() = delete;

  template <typename ...Args>
  MultiVector(Args && ... args)
    : data_( std::forward<Args>(args)... ){}

  explicit MultiVector(const wrap_t & vecobj)
    : data_(vecobj){}

  MultiVector(this_t const & other)
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

}}//end namespace rompp::core

#endif
