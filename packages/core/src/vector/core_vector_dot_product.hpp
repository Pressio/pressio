
#ifndef CORE_VECTOR_DOT_PRODUCT_HPP_
#define CORE_VECTOR_DOT_PRODUCT_HPP_

namespace core{


// template <typename wrapped_matrix_type, typename U>
// typename std::enable_if<std::is_same<U,der_t>::value >::type
// matMultiplyImpl(const matrix<wrapped_matrix_type> & matin, U& result) const
// {
//   auto & resVec = result.getNonConstRefToData();
//   auto * matPtr = matin.view();
//   assert( matPtr->cols() == this->size() );
//   //assert( result.size() == this->size() );
//   resVec = (*matPtr)*data_;
// };

// sc_t dotImpl(const der_t & b) const{
//   // // what is this?
//   // // dot product: <this,b>
//   // sc_t res = 0.0;
//   // for (size_t i=0; i<this->size(); i++)
//   //   res += data_[i]*b[i];
//   // return res;
// };

// template <typename op_t>
// void applyOpImpl(op_t op, sc_t a1,
// 		   sc_t a2, const der_t & vin){
//   // // what is this?
//   // // this = a1*this op a2*vin;
//   // for (size_t i=0; i<this->size(); i++)
//   //   data_[i] = op()( a1*data_[i], a2*vin[i] );
// }

// sc_t norm2Impl() const{
//   // sc_t result = 0;
//   // for (size_t i=0; i<this->size(); i++)
//   //   result += data_[i]*data_[i];
//   // return result;
// };

  
    
} // end namespace core
#endif

