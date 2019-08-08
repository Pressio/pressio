
#ifdef HAVE_KOKKOS
#ifndef CONTAINERS_VECTOR_CONCRETE_VECTOR_SHAREDMEM_KOKKOS_HPP_
#define CONTAINERS_VECTOR_CONCRETE_VECTOR_SHAREDMEM_KOKKOS_HPP_

#include "../../shared_base/containers_container_base.hpp"
#include "../base/containers_vector_sharedmem_base.hpp"
#include <KokkosBlas1_fill.hpp>
#include <KokkosBlas1_scal.hpp>

namespace pressio{ namespace containers{

template <typename wrapped_type>
class Vector<wrapped_type,
	     ::pressio::mpl::enable_if_t<
	       containers::meta::is_vector_kokkos<wrapped_type>::value
	       >
	     >
  : public ContainerBase< Vector<wrapped_type>, wrapped_type >,
    public VectorSharedMemBase< Vector<wrapped_type> >
{

  using this_t = Vector<wrapped_type>;
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
  Vector() = default;

  explicit Vector(const wrap_t src)
    : data_{src.label(), src.extent(0)}
  {
    Kokkos::deep_copy(data_, src);
    // std::cout << "Kokkos wrapper" << std::endl;
    // std::cout << data_.label() << std::endl;
  }

  Vector(const std::string & label, size_t e1)
    : data_{label, e1}
  {}

  Vector(const this_t & other)
    : data_{other.data_.label(), other.data_.extent(0)}
  {
    Kokkos::deep_copy(data_, other.data_);
    // std::cout << "Kokkos copy cstr" << std::endl;
    // std::cout << data_.label() << std::endl;
  }

  ~Vector(){}

public:
  // copy assign implments copy semantics not view (for time being)
  this_t & operator=(const this_t & other){
    assert(this->size() == other.size());
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

  void scaleImpl(sc_t value) {
    KokkosBlas::scal(data_, value, data_);
  }

  void setZeroImpl() {
    constexpr auto zero = ::pressio::utils::constants::zero<sc_t>();
    KokkosBlas::fill(data_, zero);
  }

  wrap_t dataCpImpl(){
    return data_;
  }

  ord_t sizeImpl() const {
    return data_.extent(0);
  }

private:
  friend ContainerBase< this_t, wrapped_type >;
  friend VectorSharedMemBase< this_t >;

private:
  wrap_t data_ = {};

};//end class

}}//end namespace pressio::containers
#endif // CONTAINERS_VECTOR_CONCRETE_VECTOR_SHAREDMEM_KOKKOS_HPP_
#endif // HAVE_TRILINOS
