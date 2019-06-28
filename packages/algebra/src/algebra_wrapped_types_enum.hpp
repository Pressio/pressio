
#ifndef ALGEBRA_WRAPPED_TYPES_ENUM_HPP_
#define ALGEBRA_WRAPPED_TYPES_ENUM_HPP_

namespace rompp{ namespace algebra{ namespace details {

/*--------------------------------------------
Wrapped library name for containers
--------------------------------------------*/
enum class WrappedPackageIdentifier{
   Undefined,
   Eigen,
   Blaze,
   Trilinos,
   Kokkos,
   Armadillo,
   Arbitrary
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
   Tpetra,
   TpetraBlock,
   Eigen,
   EigenRowStatic,
   EigenColStatic,
   EigenRowDynamic,
   EigenColDynamic,
   TeuchosSerialDense,
   BlazeStatic,
   BlazeDynamic,
   Kokkos,
   ArmadilloCol,
   ArmadilloRow,
   Arbitrary
};

/*--------------------------------------------
Identifier for wrapped matrix
--------------------------------------------*/
enum class WrappedMatrixIdentifier{
   Undefined,
   CrsEpetra,
   DenseEpetra,
   SparseEpetra,
   SparseTpetra,
   TeuchosSerialDense,
   DenseEigen, // maybe more specific, like static or dynamic
   SparseEigen,
   DenseBlaze, // maybe more specific
   SparseBlaze,
   Arbitrary
};


/*--------------------------------------------
Identifier for wrapped multivector
--------------------------------------------*/
enum class WrappedMultiVectorIdentifier{
   Undefined,
   Epetra,
   Tpetra,
   TpetraBlock,
   Eigen,
   Arbitrary
};


}}} // end namespace rompp::algebra::details
#endif
