
#ifndef CORE_VECTOR_VECTOR_TRAITS_HPP_
#define CORE_VECTOR_VECTOR_TRAITS_HPP_

#include "../core_forward_declarations.hpp"
#include "../meta/core_native_vector_meta.hpp"
#include "../meta/core_meta_detect_typedefs.hpp"
#include "../meta/core_meta_detect_operators.hpp"
#include "../core_shared_traits.hpp"

namespace rompp{
namespace core{
namespace details{


//*******************************
// Blaze dynamic vector 
//******************************* 
template <typename wrapped_type>
struct traits<Vector<wrapped_type,
		     core::meta::enable_if_t<
		       core::meta::is_dynamic_vector_blaze<
			 wrapped_type>::value
		       >
		     >
	      >
  : public containers_shared_traits<Vector<wrapped_type>,
				    wrapped_type,
				    true, false, false,
				    WrappedPackageIdentifier::Blaze,
				    true>
{

  static constexpr WrappedVectorIdentifier
  wrapped_vector_identifier = WrappedVectorIdentifier::BlazeDynamic;
    
  using scalar_t = typename wrapped_type::ElementType;
  using ordinal_t = unsigned long;

  static constexpr bool is_static = true;
  static constexpr bool is_dynamic = !is_static;
  static constexpr int rows = -1;
};
  

//*******************************
// Eigen vector 
//******************************* 
template <typename wrapped_type>
struct traits<Vector<wrapped_type,
		     typename
		     std::enable_if<
		       core::meta::is_vector_eigen<
			 wrapped_type>::value
		       >::type
		     >
	      >
  : public containers_shared_traits<Vector<wrapped_type>,
				    wrapped_type,
				    true, false, false,
				    WrappedPackageIdentifier::Eigen,
				    true>
{

  static constexpr WrappedVectorIdentifier
  wrapped_vector_identifier = WrappedVectorIdentifier::Eigen;
    
  using scalar_t = typename wrapped_type::Scalar;
  using ordinal_t = int;

  static constexpr bool is_static = (
	// if it is a row vector NON dynamic
	( wrapped_type::RowsAtCompileTime != Eigen::Dynamic &&
	  wrapped_type::ColsAtCompileTime == 1 ) ||
	// if it is a col vector NON dynamic
	( wrapped_type::RowsAtCompileTime == 1 &&
	  wrapped_type::ColsAtCompileTime != Eigen::Dynamic )
	);
  static constexpr bool is_dynamic = !is_static;
  static constexpr int rows = is_dynamic ? -1 : wrapped_type::RowsAtCompileTime;
};


//*******************************
// for epetra vector 
//******************************* 
template<typename wrapped_type>
struct traits<Vector<wrapped_type,
	  typename
	  std::enable_if<
	    core::meta::is_vector_epetra<
	      wrapped_type>::value
	    >::type
	  >
	>
  : public containers_shared_traits<Vector<wrapped_type>,
				    wrapped_type,
				    true, false, false,
			       WrappedPackageIdentifier::Trilinos,
				    false>
{

  static constexpr WrappedVectorIdentifier
  wrapped_vector_identifier = WrappedVectorIdentifier::Epetra;
  
  using scalar_t = defaultTypes::epetra_scalar_t;
  using local_ordinal_t = core::defaultTypes::epetra_lo_t;
  using global_ordinal_t = core::defaultTypes::epetra_go_t1;
  using data_map_t = Epetra_BlockMap;
  using communicator_t = Epetra_Comm;

  static constexpr bool is_dynamic = true;
  static constexpr int rows = -1;
};
  
  
//*******************************
// Kokkos vector 
//******************************* 
template <typename wrapped_type>
struct traits<Vector<wrapped_type,
	  core::meta::enable_if_t<
	    core::meta::is_vector_kokkos<
	      wrapped_type>::value
	    >
	  >
	>
  : public containers_shared_traits<Vector<wrapped_type>,
				    wrapped_type,
				    true, false, false,
				    WrappedPackageIdentifier::Kokkos,
				    true>
{

  static constexpr WrappedVectorIdentifier
  wrapped_vector_identifier = WrappedVectorIdentifier::Kokkos;

  using scalar_t = typename wrapped_type::traits::value_type;
  using ordinal_t = int;

  using  execution_space = typename wrapped_type::traits::execution_space;
  using  memory_space = typename wrapped_type::traits::memory_space;
  using  device_type = typename wrapped_type::traits::device_type;
  using  memory_traits = typename wrapped_type::traits::memory_traits;
  using  host_mirror_space = typename wrapped_type::traits::host_mirror_space;
  static constexpr bool is_static = wrapped_type::traits::rank_dynamic==0;
  static constexpr bool is_dynamic = !is_static;
  static constexpr int rows = -1; 
};


  
//*******************************
// for a std vector 
//******************************* 
template <typename wrapped_type>
struct traits<Vector<wrapped_type,
    typename
    std::enable_if<
      core::meta::is_vector_stdlib<
	wrapped_type>::value
      >::type
    >
  >
  : public containers_shared_traits<Vector<wrapped_type>,
				    wrapped_type,
				    true, false, false,
			       WrappedPackageIdentifier::CppStdLib,
				    true>
{
  
  static constexpr WrappedVectorIdentifier
  wrapped_vector_identifier = WrappedVectorIdentifier::CppStdLib;
  
  using scalar_t = typename wrapped_type::value_type;
  using ordinal_t = core::defaultTypes::local_ordinal_t;
  static constexpr bool is_dynamic = true;
  static constexpr int rows = -1;
};


  
}//end namespace details
}//end namespace core

}//end namespace rompp
#endif
