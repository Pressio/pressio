
#ifndef CORE_VECTORINTERFACE_HPP
#define CORE_VECTORINTERFACE_HPP

#include "core_ConfigDefs.hpp"
#include "forwardDeclarations.hpp"

namespace core
{
  
template<typename derived_type>
class vectorInterface
{
public:
  using sc_t = typename details::traits<derived_type>::scalar_t;
  using LO_t = typename details::traits<derived_type>::local_ordinal_t;
  using GO_t = typename details::traits<derived_type>::global_ordinal_t;
  using der_t = typename details::traits<derived_type>::derived_t;
  using wrap_t = typename details::traits<derived_type>::wrapped_t;
  using map_t = typename details::traits<derived_type>::map_t;
  using comm_t = typename details::traits<derived_type>::comm_t;
  //------------------------------------------------------------

  virtual sc_t dot (const der_t &) const = 0;
  virtual sc_t & operator [] (LO_t) = 0;
  virtual sc_t const & operator [] (LO_t) const = 0;
  virtual size_t sizeLocal() const = 0;
  virtual wrap_t const * view() const = 0;
  virtual void resize(size_t) = 0;
  virtual wrap_t & getNonConstRefToData() = 0;
  virtual sc_t norm2() const = 0;

  //-------------------
  // the following ones are defined inside impl only
  //-------------------
  // template <typename op_t>
  // void applyOp(op_t op, sc_t a1, sc_t a2, const der_t &);

  // template<class U=map_t>
  // const typename std::enable_if< !std::is_same<U,void>::value, U>::type & 
  // getMap() const{
  //   return this->underlying().getMapImpl();
  // }
  //-------------------

  
  // int replaceGlobalValue(GO_t GlobalRow, sc_t ScalarValue) = 0;
  // int replaceMyValue(LO_t MyRow, sc_t ScalarValue) = 0;
  // int abs(der_t & a) const = 0;
  // int putScalar (sc_t ScalarConstant) = 0;
  // int extractCopy(sc_t *a) const = 0;
  // int scale(sc_t ScalarValue) = 0;
  // int norm1(sc_t & Result) const = 0;
};
    
} // end namespace core

#endif



  // int replaceGlobalValue(GO_t GlobalRow, sc_t ScalarValue) override{
  //   return static_cast<der_t *>(this)->replaceGlobalValueImpl(GlobalRow, ScalarValue);
  // };
  // int replaceMyValue(LO_t MyRow, sc_t ScalarValue) override{
  //   return static_cast<der_t *>(this)->replaceMyValueImpl(MyRow, ScalarValue);
  // };
  // int abs(der_t & a) const override{
  //   return static_cast<der_t const *>(this)->absImpl(a);
  // };
  // int putScalar (sc_t ScalarConstant) override{
  //   return static_cast<der_t *>(this)->putScalarImpl(ScalarConstant);
  // };
  // int extractCopy(sc_t *a) const override{
  //   return static_cast<der_t const *>(this)->extractCopyImpl(a);
  // };
  // int scale(sc_t ScalarValue) override{
  //   return static_cast<der_t *>(this)->scaleImpl(ScalarValue);
  // };
  // int norm1(sc_t & Result) const override{
  //   return static_cast<der_t const *>(this)->norm1Impl(Result);
  // };




  
//   //! = Operator.
//   /*!
//     \param In
//     A - der_t to copy.
//     \return der_t.
//   */
//   der_t& operator = (const der_t& Source) = 0;

//   //! Vector access function.
//   const sc_t * operator [] (LO_t i) const = 0;

//   //! Returns the local vector length on the calling processor.
//   LO_t MyLength() const = 0;

//   //! Returns the global vector length.
//   GO_t GlobalLength() const = 0;
  
//   //! Replace multi-vector values with scaled values of A, \e this = ScalarA*A.
//   /*!
//     \param In
//     ScalarA - Scale value.
//     \param In
//     A - Multi-vector to copy.
//     \param Out
//     \e This - Multi-vector with values overwritten by scaled values of A.
//     \return Integer error code, set to 0 if successful.
//   */
//   int scale(sc_t ScalarA, const der_t & A) = 0;


//   //! Compute 2-norm of each vector in multi-vector.
//   /*!
//     \param Out
//     Result - Result[i] contains 2-norm of ith vector.
//     \warning Map of the \e this multivector must have unique GIDs (UniqueGIDs() must return true).
//     \return Integer error code, set to 0 if successful.
//   */
//   int norm2(sc_t & Result) const = 0;

//   //! Compute Inf-norm of each vector in multi-vector.
//   /*!
//     \param Out
//     Result - Result[i] contains Inf-norm of ith vector.
//     \return Integer error code, set to 0 if successful.
//   */
//   int normInf(sc_t * Result) const = 0;

//   //! Compute minimum value of each vector in multi-vector.
//   /*! Note that the vector contents must be already initialized for this
//       function to compute a well-defined result. The length of the
//       vector need not be greater than zero on all processors. If length is
//       greater than zero on any processor then a valid result will be computed.
//     \param Out
//     Result - Result[i] contains minimum value of ith vector.
//     \return Integer error code, set to 0 if successful.
//   */
//   int minValue(sc_t * Result) const = 0;

//   //! Compute maximum value of each vector in multi-vector.
//   /*! Note that the vector contents must be already initialized for this
//       function to compute a well-defined result. The length of the
//       vector need not be greater than zero on all processors. If length is
//       greater than zero on any processor then a valid result will be computed.
//     \param Out
//     Result - Result[i] contains maximum value of ith vector.
//     \return Integer error code, set to 0 if successful.
//   */
//   int maxValue(sc_t & Result) const = 0;

//   //! Compute mean (average) value of each vector in multi-vector.
//   /*!
//     \param Out
//     Result - Result[i] contains mean value of ith vector.
//     \warning Map of the \e this multivector must have unique GIDs (UniqueGIDs() must return true).
//     \return Integer error code, set to 0 if successful.
//   */
//   int meanValue(sc_t & Result) const = 0;

//   //! Multiply a der_t with another, element-by-element.
//   /*! This function supports diagonal matrix multiply.  A is usually a single vector
//     while B and \e this may have one or more columns.  Note that B and \e this must
//     have the same shape.  A can be one vector or have the same shape as B.  The actual
//     computation is \e this = ScalarThis * \e this + ScalarAB * B @ A where @ denotes element-wise
//     multiplication.
//   */
//   int multiply(sc_t ScalarAB, const der_t & A, const der_t & B, sc_t ScalarThis ) = 0;
// };

  


  
// template <typename T,
// 	  class sc_t, typename der_t>
// class vectorIF{
// public:
//   virtual ~vectorIF(){};
//   void dot (const der_t & ov, sc_t & result) const{
//     static_cast<der_t const *>(this)->dotImpl(ov, result);
//   };
//   T const * view() const{
//     static_cast<der_t const *>(this)->viewImpl();
//   };
// };

  
// template<class T, class sc_t>
// class vector : public vectorIF<T, sc_t, vector<T,sc_t>>{
// public:
//   T data_;
//   vector(T & vec) : data_(vec){};
//   virtual ~vector(){};
//   T const * viewImpl() const{ return &data_; };
//   void dotImpl (const vector<T, sc_t> & b, sc_t & res) const{
//     data_.dot(b.data_, res);
//   };
// };

