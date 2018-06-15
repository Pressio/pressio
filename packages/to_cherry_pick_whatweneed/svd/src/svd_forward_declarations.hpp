
#ifndef SVD_FORWARD_DECLARATIONS_HPP_
#define SVD_FORWARD_DECLARATIONS_HPP_

#include "svd_ConfigDefs.hpp"
#include "core_forward_declarations.hpp"

namespace svd {

  enum class svdKind {
    EigenJacobi,
    EigenBDCSVD,
    invalid
  };


namespace details {

  template<typename T, typename enable = void>
  struct traits; // : core::details::traits<T,enable>{};

  template<typename T>
  struct traits<const T> : traits<T>{}; //core::details::traits<const T> {};

} // end namespace details
  

  
  // forward declaration of svd class
  template <typename matrix_type,
	    svdKind which_impl = svdKind::invalid,
	    typename enable = void>
  class solver;


  // // sketch of a factory 
  // struct factory
  // {
  //   template <typename matrix_type,
  // 	      svdKind which_impl,
  // 	      typename... Args>
  //   static solver<matrix_type, which_impl> create(Args && ... args)
  //   {
  //     return solver<matrix_type, which_impl>(std::forward<Args>(args)...);
  //   };
  // };
  
  
} // end namespace 

#endif
