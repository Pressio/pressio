
#ifndef CORE_VECTORSERARB_HPP
#define CORE_VECTORSERARB_HPP


#include "core_vectorBaseImpl.hpp"
#include "core_vectorSerImpl.hpp"
#include "core_vectorMathImpl.hpp"


//***********************************
// general serial case:
// basically a user-defined vector class
//***********************************


namespace core{


  
template <typename wrapped_type>
class vector<wrapped_type,
	     typename std::enable_if<!std::is_same<wrapped_type,
						   std::vector<typename wrapped_type::scalar_type>>::value
				     >::type >
  : public vectorBaseImpl<vector<wrapped_type> >,
    public vectorSerImpl<vector<wrapped_type> >,
    public vectorMathImpl<vector<wrapped_type> >
{
public:
  using derived_t = vector<wrapped_type>;
  using sc_t = typename details::traits<derived_t>::scalar_t;
  using der_t = typename details::traits<derived_t>::derived_t;
  using wrap_t = typename details::traits<derived_t>::wrapped_t;
  using ord_t = typename details::traits<derived_t>::ordinal_t;

private:
   wrap_t data_;

public:
  vector(){};
  vector(const wrap_t & obj) : data_(obj){};
  ~vector(){};


  sc_t & operator [] (ord_t i){
    return data_[i];
  };
  sc_t const & operator [] (ord_t i) const{
    return data_[i];
  };  
  //-----------------------------------

  wrap_t const * viewImpl() const{
    return &data_;
  };
  wrap_t & getNonConstRefToDataImpl(){
    return data_;
  };
  void resizeImpl(size_t val) {
    data_.resize(val);
  };
  //-----------------------------------

  size_t sizeImpl() const {
    return data_.size();
  };
  //-----------------------------------

  sc_t dotImpl(const der_t & b) const{
    return data_.dot(b);
  };

  template <typename op_t>
  void applyOpImpl(op_t op, sc_t a1,
  		   sc_t a2, const der_t & vin){
    // what is this?
    // this = a1*this op a2*vin;
    for (size_t i=0; i<this->size(); i++)
      data_[i] = op()( a1*data_[i], a2*vin[i] );
  }
  
  sc_t norm2Impl() const{
    return data_.norm2();
  };


};
  
  
}//end namespace core
#endif
