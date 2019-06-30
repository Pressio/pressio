
#ifdef HAVE_TRILINOS
#ifndef ALGEBRA_VECTOR_CONCRETE_VECTOR_SHAREDMEM_KOKKOS_HPP_
#define ALGEBRA_VECTOR_CONCRETE_VECTOR_SHAREDMEM_KOKKOS_HPP_

#include "../../shared_base/algebra_container_base.hpp"
#include "../base/algebra_vector_sharedmem_base.hpp"

namespace rompp{ namespace algebra{

template <typename wrapped_type>
class Vector<wrapped_type,
	     ::rompp::mpl::enable_if_t<
	       algebra::meta::is_vector_kokkos<wrapped_type>::value
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

  Vector(const this_t & other)
    : data_{other.data_.label(), other.data_.extent(0)}
  {
    Kokkos::deep_copy(data_, other.data_);
    // std::cout << "Kokkos copy cstr" << std::endl;
    // std::cout << data_.label() << std::endl;
  }

  ~Vector(){}

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

private:
  friend ContainerBase< this_t, wrapped_type >;
  friend VectorSharedMemBase< this_t >;

private:
  wrap_t data_ = {};

};//end class

}}//end namespace rompp::algebra
#endif // ALGEBRA_VECTOR_CONCRETE_VECTOR_SHAREDMEM_KOKKOS_HPP_
#endif // HAVE_TRILINOS
