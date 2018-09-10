
#ifndef CORE_MULTIVECTOR_MULTIVECTOR_TRAITS_HPP_
#define CORE_MULTIVECTOR_MULTIVECTOR_TRAITS_HPP_

#include "../core_forward_declarations.hpp"
#include "../meta/core_multi_vector_meta.hpp"
#include "../core_shared_traits.hpp"

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
  : public containers_shared_traits<MultiVector<wrapped_type>,
				    wrapped_type,
				    false, false, true,
				    WrappedPackageIdentifier::Eigen,
				    false>
{

  using scalar_t = defaultTypes::epetra_scalar_t;
  using local_ordinal_t = core::defaultTypes::epetra_lo_t;
  using global_ordinal_t = core::defaultTypes::epetra_go_t1;
  using data_map_t = Epetra_BlockMap;
  using communicator_t = Epetra_Comm;

  // // a multivector can also be seen as dense matrix.
  // // But here is intended as multivector. so we set actsAsMultiVector
  // // while for dense distributed matrix, we set actsAsMultiVector = 0
  // static constexpr int actingAsMultiVector = 1;
  // static constexpr int actingAsDenseMatrix = 0; 
};

    
}//end namespace details  

  
  
namespace meta {

template <typename T, typename enable = void>
struct is_core_multi_vector_wrapper : std::false_type {};

template <typename T>
struct is_core_multi_vector_wrapper<T,
	   typename
	   std::enable_if<
	     core::details::traits<T>::is_multi_vector==1
	     >::type
	   > : std::true_type{};

#define STATIC_ASSERT_IS_CORE_MULTI_VECTOR_WRAPPER(TYPE) \
  static_assert( core::meta::is_core_multi_vector_wrapper<TYPE>::value, \
		 "THIS_IS_NOT_A_CORE_MULTI_VECTOR_WRAPPER")
  
/////////////////////////
}//end meta
/////////////////////////

  
}//end namespace core
#endif
