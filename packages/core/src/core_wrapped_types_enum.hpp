
#ifndef CORE_WRAPPED_TYPES_ENUM_HPP_
#define CORE_WRAPPED_TYPES_ENUM_HPP_

namespace rompp{
namespace core{
namespace details {  

/*--------------------------------------------
Wrapped library name for containers
--------------------------------------------*/
enum class WrappedPackageIdentifier{
   Undefined,
   Eigen,
   Trilinos,
   Kokkos,
   CppStdLib
};

/*--------------------------------------------
Identifier for wrapped vectors

Within a given package, like trilinos, we can 
have multiple types of vectors, for 
instance epetra, tpetra.
Same can be true for other packages.
--------------------------------------------*/
enum class WrappedVectorIdentifier{
   Undefined,
   Epetra,
   Eigen, // maybe more specific, like static or dynamic
   Kokkos,
   CppStdLib
};

/*--------------------------------------------
Identifier for wrapped matrix
--------------------------------------------*/
enum class WrappedMatrixIdentifier{
   Undefined,
   CrsEpetra,
   DenseEpetra,
   DenseEigen, // maybe more specific, like static or dynamic
   SparseEigen,
   CppStdLib
};


/*--------------------------------------------
Identifier for wrapped multivector
--------------------------------------------*/
enum class WrappedMultiVectorIdentifier{
   Undefined,
   Epetra,
   Eigen // maybe more specific, like static or dynamic
};

  
} // end namespace details
} // end of core namespace

}//end namespace rompp
#endif
