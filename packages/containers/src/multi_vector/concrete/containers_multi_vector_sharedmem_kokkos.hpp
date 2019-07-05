
#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_MULTI_VECTOR_CONCRETE_VECTOR_SHAREDMEM_KOKKOS_HPP_
#define CONTAINERS_MULTI_VECTOR_CONCRETE_VECTOR_SHAREDMEM_KOKKOS_HPP_

#include "../../shared_base/containers_container_base.hpp"
#include "../base/containers_multi_vector_sharedmem_base.hpp"

namespace pressio{ namespace containers{

template <typename wrapped_type>
class MultiVector<
  wrapped_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_kokkos<wrapped_type>::value
    >
  >
  : public ContainerBase< MultiVector<wrapped_type>, wrapped_type >,
    public MultiVectorSharedMemBase<MultiVector<wrapped_type>>
{

  using this_t = MultiVector<wrapped_type>;
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
  MultiVector() = default;

  explicit MultiVector(const wrap_t src)
    : data_{src.label(), src.extent(0), src.extent(1)}{
    Kokkos::deep_copy(data_, src);
  }

  MultiVector(const std::string & label, size_t e1, size_t e2)
    : data_{label, e1, e2}
  {}

  MultiVector(const this_t & other)
    : data_{other.data_.label(),
	    other.data_.extent(0),
	    other.data_.extent(1)}{
    Kokkos::deep_copy(data_, other.data_);
  }

  ~MultiVector(){}

public:
  // copy assign implments copy semantics not view (for time being)
  this_t & operator=(const this_t & other){
    assert(this->length() == other.length());
    assert(this->numVectors() == other.numVectors());
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

  ord_t numVectorsImpl() const{
    return data_.extent(1);
  }

  ord_t lengthImpl() const{
    return data_.extent(0);
  }

private:
  friend ContainerBase< this_t, wrapped_type >;
  friend MultiVectorSharedMemBase<this_t>;

private:
  wrap_t data_ = {};

};//end class

}}//end namespace pressio::containers
#endif
#endif // HAVE_TRILINOS
