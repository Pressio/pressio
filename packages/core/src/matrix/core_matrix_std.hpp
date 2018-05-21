
#ifndef CORE_MATRIX_STD_HPP_
#define CORE_MATRIX_STD_HPP_

#include "core_matrix_generic_base.hpp"
// #include "core_vector_serial_base.hpp"
// #include "core_vectorMathImpl.hpp"

// !!!!!!!!!!!!!!!!!!!!!!!!!!
// WIP
// TO FINISH 
// !!!!!!!!!!!!!!!!!!!!!!!!!!





//****************************************
// std::vector-based matrix specialization
//****************************************
namespace core{


template <typename scalar_type,
	  typename ordinal_type>
class matrix<std::vector<std::vector<scalar_type> >, scalar_type,ordinal_type>
  // : public matrixBaseImpl<matrix<std::vector<std::vector<scalar_type>>,
  // 				 scalar_type,ordinal_type> >
  //   public vectorSerialBase<vector<std::vector<scalar_type>,scalar_type,ordinal_type> >,
  //   public vectorMathImpl<vector<std::vector<scalar_type>,scalar_type,ordinal_type> >
{
public:
  using derived_t = matrix<std::vector<std::vector<scalar_type>>,scalar_type,ordinal_type>;
  using sc_t = typename details::traits<derived_t>::scalar_t;
  using der_t = typename details::traits<derived_t>::derived_t;
  using wrap_t = typename details::traits<derived_t>::wrapped_t;
  using ord_t = typename details::traits<derived_t>::ordinal_t;

private:
  std::vector<std::vector<sc_t>> data_;

public:
  matrix(){
    std::cout << "default std" << std::endl;
  }
  ~matrix(){}

  // sc_t & operator [] (ord_t i){
  //   return data_[i];
  // };
  // sc_t const & operator [] (ord_t i) const{
  //   return data_[i];
  // };  
  //-----------------------------------
  
  wrap_t const * viewImpl() const{
    return &data_;
  };
  wrap_t & getNonConstRefToDataImpl(){
    return data_;
  };
  // void resizeImpl(size_t val) {
  //   data_.resize(val);
  // };
  // //-----------------------------------

  // size_t sizeImpl() const {
  //   return data_.size();
  // };
  // //-----------------------------------
  
  // template <typename op_t>
  // void applyOpImpl(op_t op, sc_t a1,
  // 		   sc_t a2, const der_t & vin){
  //   // what is this?
  //   // this = a1*this op a2*vin;
  //   for (size_t i=0; i<this->size(); i++)
  //     data_[i] = op()( a1*data_[i], a2*vin[i] );
  // }
  
  // sc_t norm2Impl() const{
  //   sc_t result = 0;
  //   for (size_t i=0; i<this->size(); i++)
  //     result += data_[i]*data_[i];
  //   return result;
  // };


  
};
  
}//end namespace core
  
#endif

