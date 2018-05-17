
#ifndef CORE_MATRIXTRAITS_HPP
#define CORE_MATRIXTRAITS_HPP

#include "forwardDeclarations.hpp"
#include <vector>
//#include "Epetra_Vector.h"
#include <Eigen/Core>


namespace core{

namespace details{


  //*******************************
  // std::vector-based matrix wrapper 
  //******************************* 
  template <typename scalar_type,
	    typename ordinal_type>
  struct traits<matrix<std::vector<std::vector<scalar_type>>,
		       scalar_type,
		       ordinal_type> >
  {
    using scalar_t = scalar_type;
    using ordinal_t = ordinal_type;
    using wrapped_t = std::vector<std::vector<scalar_t>>;
    using derived_t = vector<wrapped_t,scalar_t,ordinal_t>;
    enum {
      isMatrix=1,
      isVector = !isMatrix,
      isSerial = 1,
      isDistributed = !isSerial
    };
  };

  //*******************************
  // eigen matrix wrapper 
  //******************************* 
  template <typename scalar_type>
  struct traits<matrix<Eigen::Matrix<scalar_type,Eigen::Dynamic,Eigen::Dynamic>,scalar_type,int>>
  {
    using scalar_t = scalar_type;
    using ordinal_t = int;
    using wrapped_t = Eigen::Matrix<scalar_type,Eigen::Dynamic,Eigen::Dynamic>;
    using derived_t = vector<wrapped_t,scalar_t>;
    enum {
      isMatrix=1,
      isVector = !isMatrix,
      isSerial = 1,
      isDistributed = !isSerial,
      isEigen = 1,
      isStatic = 0     
    };
  };

  // //****************************************
  // // eigen matrix COMPILE time size wrapper 
  // //****************************************
  // template <typename scalar_type,
  // 	    typename RowsAtCompileTime,
  // 	    typename ColsAtCompileTime>	    
  // struct traits< matrix<Eigen::Matrix<scalar_type,RowsAtCompileTime,ColsAtCompileTime>,scalar_type,int,
  // 			void,RowsAtCompileTime,ColsAtCompileTime,void,void,void
  // 			>
  // 		 >
  // {
  //   using scalar_t = scalar_type;
  //   using ordinal_t = int;
  //   using wrapped_t = Eigen::Matrix<scalar_type,RowsAtCompileTime,ColsAtCompileTime>;
  //   using derived_t = matrix<Eigen::Matrix<scalar_t,RowsAtCompileTime,ColsAtCompileTime>,
  // 			     scalar_t,int,void,RowsAtCompileTime,ColsAtCompileTime,void,void,void>;
  //   enum {
  //     isMatrix=1,
  //     isVector = !isMatrix,
  //     isSerial = 1,
  //     isDistributed = !isSerial,
  //     isEigen = 1,
  //     isStatic = 1
  //   };
  // };

  
  // //*******************************
  // // for a general distributed vector 
  // //******************************* 
  // template <typename wrapped_type,
  // 	    typename scalar_type,
  // 	    typename local_ordinal_type,
  // 	    typename global_ordinal_type,
  // 	    typename map_type,
  // 	    typename comm_type
  // 	    >
  // struct traits<vector<wrapped_type,scalar_type,local_ordinal_type,
  // 		       global_ordinal_type,map_type,comm_type>,
  // 		typename std::enable_if<!std::is_integral<scalar_type>::value &&
  // 					!std::is_same<wrapped_type,Epetra_Vector>::value && 
  // 					!std::is_void<map_type>::value &&
  // 					!std::is_void<comm_type>::value
  // 					>::type >
  // {
  //   using scalar_t = scalar_type;
  //   using local_ordinal_t = local_ordinal_type;
  //   using global_ordinal_t = global_ordinal_type;
  //   using wrapped_t = wrapped_type;
  //   using derived_t = vector<wrapped_t,scalar_t,local_ordinal_t,global_ordinal_t>;
  //   using map_t = map_type;
  //   using comm_t = comm_type;
  //   enum {
  //     isVector = 1,
  //     isSTDVector = 0,
  //     //std::is_same<wrapped_type,std::vector<scalar_t>>::value,
  //     isSerial = 0,
  //     isDistributed = 1
  //   };
  // };


  // //*******************************
  // // for epetra vector wrapper
  // //******************************* 
  // template <typename scalar_type,
  // 	    typename global_ordinal_type,
  // 	    typename map_type,
  // 	    typename comm_type
  // 	    >
  // struct traits<vector<Epetra_Vector,scalar_type,core::defaultTypes::epetra_lo_t,
  // 		       global_ordinal_type,map_type,comm_type>,
  // 		typename std::enable_if<
  // 		  std::is_same<scalar_type, core::defaultTypes::epetra_scalar_t>::value &&
  // 		  (std::is_same<global_ordinal_type,core::defaultTypes::epetra_go_t1>::value ||
  // 		   std::is_same<global_ordinal_type,core::defaultTypes::epetra_go_t2>::value) &&
  // 		  !std::is_void<map_type>::value &&
  // 		  !std::is_void<comm_type>::value
  // 		  >::type >
  // {
  //   using scalar_t = scalar_type;
  //   using local_ordinal_t = core::defaultTypes::epetra_lo_t;
  //   using global_ordinal_t = global_ordinal_type;
  //   using wrapped_t = Epetra_Vector;
  //   using map_t = map_type;
  //   using comm_t = comm_type;
  //   using derived_t = vector<wrapped_t,scalar_t,local_ordinal_t,global_ordinal_t,map_t,comm_t>;
  //   enum {
  //     isVector = 1,
  //     isSTDVector = 0,
  //     isSerial = 0,
  //     isDistributed = 1
  //   };
  // };

  
  // //*******************************
  // // for a general serial vec wrapper 
  // //******************************* 
  // template <typename wrapped_type,
  // 	    typename scalar_type,
  // 	    typename ordinal_type
  // 	    >
  // struct traits<vector<wrapped_type,scalar_type,ordinal_type>,
  // 		typename std::enable_if<
  // 		  !std::is_same<wrapped_type,std::vector<scalar_type>>::value
  // 		  >::type >
  // {
  //   using scalar_t = scalar_type;
  //   using ordinal_t = ordinal_type;
  //   using wrapped_t = wrapped_type;
  //   using derived_t = vector<wrapped_t,scalar_t,ordinal_t>;
  //   enum {
  //     isVector = 1,
  //     isSTDVector = 0,
  //     isSerial = 1,
  //     isDistributed = !isSerial
  //   };
  // };
  
  
  
}//end namespace details
}//end namespace core

#endif
