
#ifndef SVD_HACKED_EPETRA_OP_HPP_
#define SVD_HACKED_EPETRA_OP_HPP_

#include <memory>
#include "svd_solver_generic_base.hpp"
#include "CORE_ALL"

#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include <Teuchos_RCPDecl.hpp>
#include "AnasaziLOBPCGSolMgr.hpp"
#include "AnasaziSimpleLOBPCGSolMgr.hpp"
#include "Epetra_InvOperator.h"


namespace svd{

template<typename matrix_type,
	 template <typename...> class lsv_type,
	 template <typename...> class rsv_type>
class Solver<matrix_type, lsv_type, rsv_type,
	     typename
	     std::enable_if<
	       core::meta::is_core_matrix_wrapper<matrix_type>::value==true &&
	       core::details::traits<matrix_type>::isEpetra==1 &&
	       core::details::traits<matrix_type>::isSparse==1
	       >::type
	     >
    : public SolverBase< Solver<matrix_type, lsv_type, rsv_type> >
{

private:
  using MV = Epetra_MultiVector;
  using OP = Epetra_Operator;
  using sc_t = typename core::details::traits<matrix_type>::scalar_t;
  using lsv_t = lsv_type<MV>;
  using rsv_t = rsv_type<MV>;
  
  const sc_t tol = 1e-10;
  const int maxIters = 500;
  const int blockSize = 2;
  
public:
  Solver() = default;
  ~Solver() = default;
  
private:

  void computeImpl(matrix_type & A, int nlsv, int nrsv){

    using Teuchos::RCP;
    
    // check that number of vectors requested is admissible
    auto nrA = A.globalRows();
    auto ncA = A.globalCols();
    assert( nlsv <= nrA );
    assert( nrsv <= ncA );

    assert(A.isFillingCompleted());
    
    //    auto ASWT = core::transpose(A);
    if (nlsv > 0){
      //***********************************************
      // compute left-sing vectors as solution of
      // A*A^T U = lambda U
      //***********************************************
      //auto AAT = core::matrixMatrixProduct(A, ASWT, false, false, true);
      auto AAT = core::matrixMatrixProduct(A, A, false, true, true);
      assert( AAT.isFillingCompleted() );

      // wrap the OP into RCP: prefer rcpFromRef
      RCP<Epetra_CrsMatrix> Arcp = Teuchos::rcpFromRef(*AAT.data() );
      // //Teuchos::RCP<Epetra_CrsMatrix> Arcp =
      // // Teuchos::rcp(new Epetra_CrsMatrix(*AAT.data()));

      ///solve eigv problem
      RCP<MV> U = this->solveEigProblem(Arcp, nlsv);
      lsv_ = std::make_shared<lsv_t>(*U);
    }

    if (nrsv > 0 ){
      //***********************************************
      // compute right-sing vectors as solution of
      // A^T * A V = lambda V
      //***********************************************
      //      auto ATA = core::matrixMatrixProduct(ASWT, A, false, false, true);
      auto ATA = core::matrixMatrixProduct(A, A, true, false, true);

      /// need to wrap the OP into RCP: prefer rcpFromRef
      RCP<Epetra_CrsMatrix> Arcp2 = Teuchos::rcpFromRef( *ATA.data() );

      // solve eigv problem
      Teuchos::RCP<MV> V = this->solveEigProblem(Arcp2, nrsv);
      rsv_ = std::make_shared<rsv_t>(*V);
    }

  }//end method
  //-----------------------------------------------------

  
  const lsv_t & cRefLeftSingularVectorsImpl() const {
    return *lsv_;
  }//end method
  //-----------------------------------------------------

  
  const rsv_t & cRefRightSingularVectorsImpl() const {
    return *rsv_;
  }//end method
  //-----------------------------------------------------


  
private:
  
  auto solveEigProblem(Teuchos::RCP<Epetra_CrsMatrix> Arcp, int nev)
  {
    auto & dommap = Arcp->DomainMap();
    Teuchos::RCP<MV> ivec = Teuchos::rcp(new MV(dommap, blockSize));
    ivec->Random();

    //-----------------------------------------------------------------
    // Create the eigenproblem. 
    Teuchos::RCP<Anasazi::BasicEigenproblem<sc_t, MV, OP> > problem =
      Teuchos::rcp(new Anasazi::BasicEigenproblem<sc_t,MV,OP>(Arcp,ivec));
    problem->setNEV(nev);
    problem->setHermitian(true);
    ///  problem->setPrec (PrecOp);
    ///finishing passing it information.
    assert(problem->setProblem());
  
    //-----------------------------------------------------------------
    // Create the solver manager and solve
    Teuchos::ParameterList PL;
    PL.set ("Which", "LM"); // ordering of eigenvalues, lowest to largest
    PL.set ("Block Size", blockSize);
    PL.set ("Maximum Iterations", maxIters);
    PL.set ("Convergence Tolerance", tol);
    const int verbosity = Anasazi::Errors + \
      /*Anasazi::Warnings +*/ Anasazi::FinalSummary;
    PL.set ("Verbosity", verbosity);

    // actual solver object
    Anasazi::SimpleLOBPCGSolMgr<sc_t, MV, OP> solverMan(problem, PL);
    auto returnCode = solverMan.solve();
    assert(returnCode == Anasazi::Converged );

    //-----------------------------------------------------------------
    // get the solution
    const Anasazi::Eigensolution<sc_t,MV> & sol = problem->getSolution();
    //const int numCompEV = sol.numVecs;
    std::vector<Anasazi::Value<sc_t> > evals = sol.Evals;
    
    std::transform(evals.cbegin(), evals.cend(),
    		   std::back_inserter(eigVals_),
    		   [](const Anasazi::Value<sc_t> & value) {
    		     return std::sqrt(value.realpart);
    		   });
    for (auto it : eigVals_)
      std::cout << it << " "  << std::endl;  
	 
    return sol.Evecs;
  }// end method
  //-----------------------------------

 
private:
  friend SolverBase< Solver<matrix_type, lsv_type, rsv_type> >;

  std::shared_ptr<lsv_t> lsv_;
  std::shared_ptr<rsv_t> rsv_;
  std::vector<sc_t> eigVals_;
  
};//end class

  
}//end namespace svd
#endif 
