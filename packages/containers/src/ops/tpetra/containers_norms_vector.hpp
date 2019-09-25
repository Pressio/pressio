
#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_SRC_OPS_TPETRA_NORMS_HPP_
#define CONTAINERS_SRC_OPS_TPETRA_NORMS_HPP_

#include "../containers_ops_meta.hpp"
#include "../../vector/containers_vector_meta.hpp"

namespace pressio{ namespace containers{ namespace ops{

template <
  typename vec_type,
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

}}}//end namespace pressio::containers::ops
#endif
#endif
