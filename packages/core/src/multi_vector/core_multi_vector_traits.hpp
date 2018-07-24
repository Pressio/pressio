
#ifndef CORE_MULTIVECTOR_MULTIVECTOR_TRAITS_HPP_
#define CORE_MULTIVECTOR_MULTIVECTOR_TRAITS_HPP_

#include "../core_forward_declarations.hpp"
#include "../meta/core_multi_vector_meta.hpp"

namespace core{
namespace details{

//*******************************
// for epetra multivector 
//******************************* 
template<typename wrapped_type>
struct traits<MultiVector<wrapped_type,
      typename std::enable_if<
       meta::is_multi_vector_epetra<wrapped_type
      >::value>::type>
     >
{

  using scalar_t = defaultTypes::epetra_scalar_t;
  using local_ordinal_t = core::defaultTypes::epetra_lo_t;
  using global_ordinal_t = core::defaultTypes::epetra_go_t1;
  using wrapped_t = Epetra_MultiVector;
  using data_map_t = Epetra_BlockMap;
  using communicator_t = Epetra_Comm;
  using derived_t = MultiVector<wrapped_t>;
  static constexpr int isVector = 0;
  static constexpr int isMultiVector = 1;
  static constexpr int isDistributed = 1;
  static constexpr int isEpetra = 1;
  static constexpr int isSTDVector = 0;
  static constexpr int isSerial = 0;
  static constexpr int isEigen = 0;
  // make these void just to be clear they are not usable
  using ordinal_t = void;
};

  
  
/////////////////////////
}//end namespace details  
/////////////////////////

  
namespace meta {

template <typename T, typename enable = void>
struct is_coreMultiVector : std::false_type {};

template <typename T>
struct is_coreMultiVector<T,
	   typename
	   std::enable_if<
	     core::details::traits<T>::isMultiVector==1
	     >::type
	   > : std::true_type{};

// #define STATIC_ASSERT_IS_CORE_MULTI_VECTOR_WRAPPER(TYPE) \
//   static_assert( core::meta::is_coreMultiVector<TYPE>::value, \
// 		 "THIS_IS_NOT_A_CORE_MULTI_VECTOR_WRAPPER")

/////////////////////////
}//end meta
/////////////////////////

  
}//end namespace core
#endif
