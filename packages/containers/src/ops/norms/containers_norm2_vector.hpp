
#ifndef CONTAINERS_CONTAINER_OPS_NORMS_NORM2_VECTOR_HPP_
#define CONTAINERS_CONTAINER_OPS_NORMS_NORM2_VECTOR_HPP_

#include "../containers_ops_meta.hpp"
#include "../../vector/containers_vector_meta.hpp"
#ifdef HAVE_KOKKOS
#include "KokkosBlas1_nrm2.hpp"
#endif

namespace pressio{ namespace containers{ namespace ops{

#ifdef HAVE_TRILINOS
//  epetra vector wrapper
template <typename vec_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_epetra<vec_type>::value
    > * = nullptr
  >
auto norm2(const vec_type & a)
  -> typename details::traits<vec_type>::scalar_t
{
  using sc_t = typename details::traits<vec_type>::scalar_t;
  sc_t result = 0.0;
  a.data()->Norm2(&result);
  return result;
}

//  tpetra vector wrapper
template <typename vec_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_tpetra<vec_type>::value &&
    std::is_same<typename details::traits<vec_type>::scalar_t,
		 typename details::traits<vec_type>::mag_t>::value
    > * = nullptr
  >
auto norm2(const vec_type & a)
  -> typename details::traits<vec_type>::mag_t
{
  return a.data()->norm2();
}

//  block tpetra vector wrapper
template <typename vec_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_tpetra_block<vec_type>::value &&
    std::is_same<typename details::traits<vec_type>::scalar_t,
		 typename details::traits<vec_type>::mag_t>::value
    > * = nullptr
  >
auto norm2(const vec_type & a)
  -> typename details::traits<vec_type>::mag_t
{
  /* workaround the non-constness of getVectorView,
   * which is supposed to be const but it is not */
  using mv_t = Tpetra::Experimental::BlockVector<>;
  return const_cast<mv_t*>(a.data())->getVectorView().norm2();
}

//  teuchos serial dense vector wrapper
template <typename vec_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_dense_vector_wrapper_teuchos<vec_type>::value
    > * = nullptr
  >
auto norm2(const vec_type & a)
  -> typename details::traits<vec_type>::scalar_t
{
  using sc_t = typename details::traits<vec_type>::scalar_t;
  sc_t result = 0.0;
  for (decltype(a.size()) i=0; i<a.size(); i++)
    result += a[i]*a[i];
  return std::sqrt(result);
}
#endif


#ifdef HAVE_KOKKOS
//  kokkos vector wrapper
template <typename vec_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_kokkos<vec_type>::value
    > * = nullptr
  >
auto norm2(const vec_type & a)
  -> typename details::traits<vec_type>::scalar_t
{
  return KokkosBlas::nrm2(*a.data());
}
#endif
//--------------------------------------------------------


//  eigen vector wrapper
template <typename vec_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_eigen<vec_type>::value
    > * = nullptr
  >
auto norm2(const vec_type & a)
  -> typename details::traits<vec_type>::scalar_t
{
  using sc_t = typename details::traits<vec_type>::scalar_t;
  sc_t result = 0.0;
  for (decltype(a.size()) i=0; i<a.size(); i++)
    result += a[i]*a[i];
  return std::sqrt(result);
}
//--------------------------------------------------------

#ifdef HAVE_PYBIND11
// pybind11::array
template <typename vec_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_array_pybind11<vec_type>::value
    > * = nullptr
  >
auto norm2(const vec_type & a)
  -> typename vec_type::value_type
{
  using sc_t = typename vec_type::value_type;
  sc_t result = 0.0;

  // make sure this is a vector
  if (a.ndim() > 1){
    throw std::runtime_error("a.ndims()!=1, this norm op is for pybind11 vectors");
  }

  for (decltype(a.size()) i=0; i<a.size(); i++)
    result += a.at(i)*a.at(i);
  return std::sqrt(result);
}
#endif
//--------------------------------------------------------

}}}//end namespace pressio::containers::ops
#endif
