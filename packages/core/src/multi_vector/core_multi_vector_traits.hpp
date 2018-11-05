
#ifndef CORE_MULTIVECTOR_MULTIVECTOR_TRAITS_HPP_
#define CORE_MULTIVECTOR_MULTIVECTOR_TRAITS_HPP_

#include "../core_forward_declarations.hpp"
#include "../meta/core_native_multi_vector_meta.hpp"
#include "../core_shared_traits.hpp"

namespace rompp{ namespace core{ namespace details{

#ifdef HAVE_TRILINOS 
//*******************************
// for epetra multivector 
//******************************* 
template<typename wrapped_type>
struct traits<MultiVector<wrapped_type,
      typename std::enable_if<
       meta::is_multi_vector_epetra<wrapped_type
      >::value>::type>
     >
  : public containers_shared_traits<MultiVector<wrapped_type>,
				    wrapped_type,
				    false, false, true,
				    WrappedPackageIdentifier::Trilinos,
				    false>
{
  static constexpr WrappedMultiVectorIdentifier
  wrapped_multi_vector_identifier = WrappedMultiVectorIdentifier::Epetra;

  using scalar_t = defaultTypes::epetra_scalar_t;
  using local_ordinal_t = core::defaultTypes::epetra_lo_t;
  using global_ordinal_t = core::defaultTypes::epetra_go_t1;
  using data_map_t = Epetra_BlockMap;
  using communicator_t = Epetra_Comm;
};
#endif


#ifdef HAVE_TRILINOS 
//*******************************
// for tpetra multivector 
//******************************* 
template<typename wrapped_type>
struct traits<MultiVector<wrapped_type,
      typename std::enable_if<
       meta::is_multi_vector_tpetra<wrapped_type
      >::value>::type>
     >
  : public containers_shared_traits<MultiVector<wrapped_type>,
				    wrapped_type,
				    false, false, true,
				    WrappedPackageIdentifier::Trilinos,
				    false>
{
  static constexpr WrappedMultiVectorIdentifier
  wrapped_multi_vector_identifier = WrappedMultiVectorIdentifier::Tpetra;

  using scalar_t = typename wrapped_type::impl_scalar_type;
  using local_ordinal_t = typename wrapped_type::local_ordinal_type;
  using global_ordinal_t = typename wrapped_type::global_ordinal_type;
  using data_map_t = typename wrapped_type::map_type;
  using device_t = typename wrapped_type::device_type;
  using node_t = typename wrapped_type::node_type;
  using dual_view_t = typename wrapped_type::dual_view_type;
  using dot_t = typename wrapped_type::dot_type;
  using mag_t = typename wrapped_type::mag_type;
  using communicator_t = decltype(std::declval<data_map_t>().getComm());
};
#endif
      

//*******************************
// for eigen multivector 
//******************************* 
template<typename wrapped_type>
struct traits<MultiVector<wrapped_type,
      typename std::enable_if<
       meta::is_multi_vector_eigen_dynamic<wrapped_type
      >::value>::type>
     >
  : public containers_shared_traits<MultiVector<wrapped_type>,
            wrapped_type,
            false, false, true,
            WrappedPackageIdentifier::Eigen,
            true>
{
  static constexpr WrappedMultiVectorIdentifier
  wrapped_multi_vector_identifier = WrappedMultiVectorIdentifier::Eigen;

  using scalar_t = typename wrapped_type::Scalar;
  using ordinal_t = int;

  static constexpr bool is_static =
    ( wrapped_type::RowsAtCompileTime != Eigen::Dynamic &&
      wrapped_type::ColsAtCompileTime != Eigen::Dynamic );
  static constexpr bool is_dynamic = !is_static;
};
    
}}}//end namespace rompp::core::details
#endif
