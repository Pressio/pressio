
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
   Blaze,
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
   BlazeStatic,
   BlazeDynamic,
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
   SparseEpetra,
   DenseEigen, // maybe more specific, like static or dynamic
   SparseEigen,
   DenseBlaze, // maybe more specific
   SparseBlaze,
   CppStdLib
};


/*--------------------------------------------
Identifier for wrapped multivector
--------------------------------------------*/
enum class WrappedMultiVectorIdentifier{
   Undefined,
   Epetra
};

  
} // end namespace details
} // end of core namespace

}//end namespace rompp
#endif
