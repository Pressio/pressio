
#ifndef CORE_MULTIVECTOR_MULTIVECTOR_TRAITS_HPP_
#define CORE_MULTIVECTOR_MULTIVECTOR_TRAITS_HPP_

#include "../core_forward_declarations.hpp"
#include "../meta/core_native_multi_vector_meta.hpp"
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
    
}//end namespace details
}//end namespace core
#endif
