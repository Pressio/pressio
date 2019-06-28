
#ifdef HAVE_TRILINOS
#ifndef ALGEBRA_VECTOR_CONCRETE_VECTOR_SHAREDMEM_KOKKOS_HPP_
#define ALGEBRA_VECTOR_CONCRETE_VECTOR_SHAREDMEM_KOKKOS_HPP_

#include "../../shared_base/algebra_container_base.hpp"
//#include "../../shared_base/algebra_container_resizable_base.hpp"
//#include "../../shared_base/algebra_container_nonresizable_base.hpp"
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

  // Views have "view semantics."  This means that they behave like
  // pointers, not like std::vector.  Their copy constructor and
  // operator= only do shallow copies.  Thus, you can pass View
  // objects around by "value"; they won't do a deep copy unless you
  // explicitly ask for a deep copy.

public:
  Vector() = default;

  explicit Vector(const wrap_t src) // by value
    : data_(src){}

  ~Vector(){}

private:

  // wrap_t const * dataImpl() const{
  //   return &data_;
  // }
  // wrap_t * dataImpl(){
  //   return &data_;
  // }

  //
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
