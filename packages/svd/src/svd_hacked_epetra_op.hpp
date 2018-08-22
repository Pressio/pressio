
#ifndef SVD_HACKED_EPETRA_OP_HPP_
#define SVD_HACKED_EPETRA_OP_HPP_

#include "svd_ConfigDefs.hpp"
#include "svd_solver_generic_base.hpp"

namespace svd{

template<typename matrix_type,
	 typename precond_type>
class Solver<matrix_type,
	     precond_type,
	     typename
	     std::enable_if<
	       core::meta::is_core_matrix_wrapper<matrix_type>::value==true &&
	       core::details::traits<matrix_type>::isEpetra==1 
	       >::type
	     >
  : public SolverBase< Solver<matrix_type> >
{

private:
  using Teuchos::RCP;
  using Teuchos::rcp;    
  using MV = Epetra_MultiVector;
  using OP = Epetra_Operator;

public:
  Solver() = default;
  ~Solver() = default;
  
private:

  template <typename U = matrix_type,
  	    typename
  	    std::enable_if<
  	      core::details::traits<U>::isSparse==1
  	      >::type * = nullptr
  	    >
  void computeImpl(const U & A, bool isHermitian, int nev)
  {
    using sc_t = typename core::details::traits<matrix_type>::scalar_t;
    using Teuchos::RCP;
    using Teuchos::rcp;    
    using MV = Epetra_MultiVector;
    using OP = Epetra_Operator;
    
    //-----------------------------------------------  
    // create data structures
    // need to wrap the OP into RCP
    //RCP<const Epetra_CrsMatrix> Arcp = Teuchos::rcpFromRef( *A.data() );
    RCP<Epetra_CrsMatrix> Arcp = rcp( new Epetra_CrsMatrix(*A.data()) );
    auto & rowmap = A.getRowDataMap();
    RCP<MV> ivec = rcp(new MV(rowmap, blockSize));
    ivec->Random();


    
  }
  //-----------------------------------------------------


  
private:
  friend SolverBase< Solver<matrix_type, precond_type> >;

  RCP<MV> evecs_;
  
};//end class

}//end namespace svd
#endif 
