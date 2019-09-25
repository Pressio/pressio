
#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_SRC_OPS_EPETRA_MULTI_VECTOR_PROD_VECTOR_HPP_
#define CONTAINERS_SRC_OPS_EPETRA_MULTI_VECTOR_PROD_VECTOR_HPP_

#include "../containers_ops_meta.hpp"
#include "../../multi_vector/containers_multi_vector_meta.hpp"

namespace pressio{ namespace containers{ namespace ops{

/*
 * multi_vector prod vector
 */


//-------------------------------------------------------//
//  EPETRA multivector with eigen vector
//-------------------------------------------------------//

// the result type is an Epetra wrapper and object is passed in
template <
  typename mvec_type,
  typename vec_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_epetra<mvec_type>::value and
    containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value and
    containers::meta::is_vector_wrapper_eigen<vec_type>::value
    > * = nullptr
  >
void product(const mvec_type & mvA,
	     const vec_type & vecB,
	     containers::Vector<Epetra_Vector> & C){

  //zero out result
  C.setZero();
  // how many vectors are in mvA
  const auto numVecs = mvA.globalNumVectors();
  // size of vecB
  assert(size_t(numVecs) == size_t(vecB.size()));
  // the data map of the multivector
  const auto mvMap = mvA.getDataMap();
  // my number of rows
  const auto myNrows = mvMap.NumMyElements();

  // loop
  for (size_t i=0; i<(size_t)myNrows; i++){
    for (size_t j=0; j<(size_t)numVecs; j++){
      C[i] += mvA(i,j) * vecB[j];
    }
  }
}

// return result Epetra wrapper
template <
  typename mvec_type,
  typename vec_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_epetra<mvec_type>::value and
    containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value and
    containers::meta::is_vector_wrapper_eigen<vec_type>::value
  > * = nullptr
 >
containers::Vector<Epetra_Vector>
product(const mvec_type & mvA, const vec_type & vecB) {

  // here, mvA is distrubted, but vecB is NOT.
  // we interpret this as a linear combination of vectors

  // the data map of the multivector
  const auto mvMap = mvA.getDataMap();
  // result is an Epetra Vector with same distribution of mvA
  using res_t = containers::Vector<Epetra_Vector>;
  res_t c(mvMap);
  product(mvA, vecB, c);
  return c;
}


//-------------------------------------------------------//
//  EPETRA multivector with teuchos vector
//-------------------------------------------------------//

template <
  typename mvec_type,
  typename vec_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_epetra<mvec_type>::value and
    containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value and
    containers::meta::is_dense_vector_wrapper_teuchos<vec_type>::value
    > * = nullptr
  >
void product(const mvec_type & mvA,
	     const vec_type & vecB,
	     containers::Vector<Epetra_Vector> & C){

  //zero out result
  C.setZero();
  // how many vectors are in mvA
  const auto numVecs = mvA.globalNumVectors();
  // size of vecB
  assert(size_t(numVecs) == size_t(vecB.size()));
  // the data map of the multivector
  const auto mvMap = mvA.getDataMap();
  // my number of rows
  const auto myNrows = mvMap.NumMyElements();

  // loop
  for (decltype(myNrows) i=0; i<myNrows; i++){
    for (decltype(numVecs) j=0; j<numVecs; j++){
      C[i] += mvA(i,j) * vecB[j];
    }
  }
}

// result is returned
template <
  typename mvec_type,
  typename vec_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_epetra<mvec_type>::value and
    containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value and
    containers::meta::is_dense_vector_wrapper_teuchos<vec_type>::value
  > * = nullptr
 >
containers::Vector<Epetra_Vector>
product(const mvec_type & mvA, const vec_type & vecB) {

  // here, mvA is distrubted, but vecB is NOT.
  // we interpret this as a linear combination of vectors

  // the data map of the multivector
  const auto mvMap = mvA.getDataMap();
  // result is an Epetra Vector with same distribution of mvA
  using res_t = containers::Vector<Epetra_Vector>;
  res_t c(mvMap);
  product(mvA, vecB, c);
  return c;
}





//-------------------------------------------------------//
//  EPETRA multivector with armadillo vector
//-------------------------------------------------------//
// //-----------------------------------------------------
// // we pass the result object
// template <
//   typename mvec_type,
//   typename vec_type,
//   ::pressio::mpl::enable_if_t<
//     containers::meta::is_multi_vector_wrapper_epetra<mvec_type>::value and
//     containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value and
//     (containers::meta::is_column_vector_wrapper_armadillo<vec_type>::value or
//      containers::meta::is_row_vector_wrapper_armadillo<vec_type>::value)
//     > * = nullptr
//   >
// void product(const mvec_type & mvA,
// 	     const vec_type & vecB,
// 	     containers::Vector<Epetra_Vector> & C){

//   //zero out result
//   C.setZero();
//   // how many vectors are in mvA
//   const auto numVecs = mvA.globalNumVectors();
//   // size of vecB
//   const size_t vecBLen = vecB.size();
//   assert(size_t(numVecs) == vecBLen);
//   // the data map of the multivector
//   const auto mvMap = mvA.getDataMap();
//   // my number of rows
//   const auto myNrows = mvMap.NumMyElements();

//   // loop
//   for (decltype(myNrows) i=0; i<myNrows; i++){
//     for (decltype(numVecs) j=0; j<numVecs; j++){
//       C[i] += mvA(i,j) * vecB[j];
//     }
//   }
// }
// //-------------------------------------------------------

// // result is returned
// template <typename mvec_type,
// 	  typename vec_type,
//   ::pressio::mpl::enable_if_t<
//     containers::meta::is_multi_vector_wrapper_epetra<mvec_type>::value and
//     containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value and
//     (containers::meta::is_column_vector_wrapper_armadillo<vec_type>::value or
//      containers::meta::is_row_vector_wrapper_armadillo<vec_type>::value)
//   > * = nullptr
//  >
// containers::Vector<Epetra_Vector>
// product(const mvec_type & mvA, const vec_type & vecB) {

//   // here, mvA is distrubted, but vecB is NOT.
//   // we interpret this as a linear combination of vectors

//   // the data map of the multivector
//   const auto mvMap = mvA.getDataMap();
//   // result is an Epetra Vector with same distribution of mvA
//   using res_t = containers::Vector<Epetra_Vector>;
//   res_t c(mvMap);
//   product(mvA, vecB, c);
//   return c;
// }


}}}//end namespace pressio::containers::ops
#endif
#endif
