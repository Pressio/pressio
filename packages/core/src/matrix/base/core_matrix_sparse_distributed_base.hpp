
#ifndef CORE_MATRIX_SPARSE_DISTRIBUTED_BASE_HPP_
#define CORE_MATRIX_SPARSE_DISTRIBUTED_BASE_HPP_

#include "../core_matrix_traits.hpp"

namespace core{
    
template<typename derived_type>
class matrixSparseDistributedBase
{

private:
  using traits_t = details::traits<derived_type>;
  using der_t = typename traits_t::derived_t;
  using wrap_t = typename traits_t::wrapped_t;

  using sc_t = typename traits_t::scalar_t;
  using LO_t = typename traits_t::local_ordinal_t;
  using GO_t = typename traits_t::global_ordinal_t;
  using row_map_t = typename traits_t::row_map_t;
  using col_map_t = typename traits_t::col_map_t;
  using comm_t =  typename traits_t::communicator_t;
  
public:
  void insertGlobalValues(GO_t targetRow,
			  LO_t numEntries,
			  const sc_t * values,
			  const GO_t * indices){
    this->underlying().insertGlobalValuesImpl(targetRow,
					      numEntries,
					      values,
					      indices);
  }

private:  
  friend der_t;
  matrixSparseDistributedBase() = default;
  ~matrixSparseDistributedBase() = default; 

private:  
  der_t & underlying(){
    return static_cast<der_t &>(*this);
  };

  der_t const& underlying() const{
    return static_cast<der_t const&>(*this);
  };

};//end class
} // end namespace core
#endif
