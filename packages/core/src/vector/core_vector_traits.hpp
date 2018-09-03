
#ifndef CORE_VECTOR_VECTOR_TRAITS_HPP_
#define CORE_VECTOR_VECTOR_TRAITS_HPP_

#include "../core_forward_declarations.hpp"
#include "../meta/core_vector_meta.hpp"
#include "../meta/core_meta_detect_typedefs.hpp"
#include "../meta/core_meta_detect_operators.hpp"

namespace core{
namespace details{

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
{

  using scalar_t = typename wrapped_type::Scalar;
  using ordinal_t = int;
  using wrapped_t = wrapped_type;
  using derived_t = Vector<wrapped_t>;
  static constexpr int isVector = 1;
  static constexpr int isEigen = 1;
  static constexpr int isSharedMem = 1;
  static constexpr int isSTDVector = 0;
  static constexpr int isDistributed = 0;
  static constexpr int isStatic = (
	// if it is a row vector NON dynamic
	( wrapped_t::RowsAtCompileTime != Eigen::Dynamic &&
	  wrapped_t::ColsAtCompileTime == 1 ) ||
	// if it is a col vector NON dynamic
	( wrapped_t::RowsAtCompileTime == 1 &&
	  wrapped_t::ColsAtCompileTime != Eigen::Dynamic )
	);
  // make these void just to be clear they are not usable
  using local_ordinal_t = void;
  using global_ordinal_t = void;
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
{

  using scalar_t = typename wrapped_type::value_type;
  using ordinal_t = core::defaultTypes::local_ordinal_t;
  using wrapped_t = wrapped_type;
  using derived_t = Vector<wrapped_t>;
  static constexpr int isVector = 1;
  static constexpr int isSTDVector = 1;
  static constexpr int isSharedMem = 1;
  static constexpr int isDistributed = 0;
  static constexpr int isEigen = 0;
  // make these void just to be clear they are not usable
  using local_ordinal_t = void;
  using global_ordinal_t = void;
};


//*******************************
// user-defined sharedMem vector 
//******************************* 
template <typename wrapped_type>
struct traits<Vector<wrapped_type,
      typename
      std::enable_if<
           core::meta::is_detected<
	     core::meta::has_scalar_typedef, wrapped_type
	     >::value &&
	  core::meta::has_ordinal_typedef<
	     wrapped_type>::value &&
	  !core::meta::is_vector_stdlib<
	     wrapped_type>::value &&
	  !core::meta::is_vector_eigen<
	     wrapped_type>::value
	  >::type
	>
     >
{

  using scalar_t = typename wrapped_type::scalar_type;
  using ordinal_t = typename wrapped_type::ordinal_type;
  using wrapped_t = wrapped_type;
  using derived_t = Vector<wrapped_type>;
  static constexpr int isVector = 1;
  static constexpr int isSharedMem = 1;
  static constexpr int isDistributed = 0;
  static constexpr int isEigen = 0;
  static constexpr int isSTDVector = 0;
  // make these void just to be clear they are not usable
  using local_ordinal_t = void;
  using global_ordinal_t = void;

};


//*******************************
// for epetra vector 
//******************************* 
template<typename wrapped_type>
struct traits<Vector<wrapped_type,
      typename
      std::enable_if<
      core::meta::is_vector_epetra<wrapped_type
      >::value>::type>
      >
{
  
  using scalar_t = defaultTypes::epetra_scalar_t;
  using local_ordinal_t = core::defaultTypes::epetra_lo_t;
  using global_ordinal_t = core::defaultTypes::epetra_go_t1;
  using wrapped_t = Epetra_Vector;
  using data_map_t = Epetra_BlockMap;
  using communicator_t = Epetra_Comm;
  using derived_t = Vector<wrapped_t>;
  static constexpr int isVector = 1;
  static constexpr int isMultiVector = 0;
  static constexpr int isDistributed = 1;
  static constexpr int isEpetra = 1;
  static constexpr int isSTDVector = 0;
  static constexpr int isSharedMem = 0;
  static constexpr int isEigen = 0;
  // make these void just to be clear they are not usable
  using ordinal_t = void;

};

  
/////////////////////////
}//end namespace details  
/////////////////////////

  
namespace meta {

template <typename T, typename enable = void>
struct is_core_vector : std::false_type {};

template <typename T>
struct is_core_vector<T,
	   typename
	   std::enable_if<
	     core::details::traits<T>::isVector==1
	     >::type
	   > : std::true_type{};

#define STATIC_ASSERT_IS_CORE_VECTOR_WRAPPER(TYPE) \
  static_assert( core::meta::is_core_vector<TYPE>::value, \
		 "THIS_IS_NOT_A_CORE_VECTOR_WRAPPER")

/////////////////////////
}//end meta
/////////////////////////

  
}//end namespace core
#endif
