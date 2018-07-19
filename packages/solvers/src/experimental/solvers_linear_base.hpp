
#ifndef SOLVERS_EXPERIMENTAL_LINEAR_BASE_HPP_
#define SOLVERS_EXPERIMENTAL_LINEAR_BASE_HPP_

#include "solvers_ConfigDefs.hpp"
#include "solvers_forward_declarations.hpp"


namespace solvers{
namespace experimental{

  
template<typename derived_type,
	 typename matrix_type,
	 typename rhs_type,
	 typename result_type
	 >
class linearBase{
public:  
  void solve(const matrix_type & A,
	     const rhs_type & b,
	     result_type & x){
    this->underlying().solveImpl(A,b,x);
  }

private:
  friend derived_type;
  linearBase() = default;
  ~linearBase() = default;

  derived_type & underlying(){
    return static_cast<derived_type& >( *this );
  }
  derived_type const & underlying() const{
    return static_cast<derived_type const & >( *this );
  } 
  
};//end class
}//end namespace experimental
}//end namespace solvers
#endif
