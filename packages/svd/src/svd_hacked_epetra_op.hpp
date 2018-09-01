
#ifndef SVD_HACKED_EPETRA_OP_HPP_
#define SVD_HACKED_EPETRA_OP_HPP_

//#include <memory>
#include "svd_solver_generic_base.hpp"
#include "CORE_ALL"

#include "AnasaziConfigDefs.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBlockDavidsonSolMgr.hpp"
#include "AnasaziSimpleLOBPCGSolMgr.hpp"
#include "AnasaziLOBPCGSolMgr.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_InvOperator.h"
#include "Ifpack.h"

namespace svd{

template<typename matrix_type,
	 template <typename...> class lsv_type,
	 template <typename...> class rsv_type,
	 typename sval_type>
class Solver<matrix_type, lsv_type, rsv_type, sval_type,
	     typename
	     std::enable_if<
	       core::meta::is_core_matrix_wrapper<matrix_type>::value==true &&
	       core::details::traits<matrix_type>::isEpetra==1 &&
	       core::details::traits<matrix_type>::isSparse==1
	       >::type
	     >
  : public SolverBase< Solver<matrix_type, lsv_type,
			      rsv_type, sval_type> >
{

private:
  using MV = Epetra_MultiVector;
  using OP = Epetra_Operator;
  using sc_t = typename core::details::traits<matrix_type>::scalar_t;
  using lsv_t = lsv_type<MV>;
  using rsv_t = rsv_type<MV>;
  using sval_t = sval_type;
  using eigprob_t = Teuchos::RCP<Anasazi::BasicEigenproblem<sc_t, MV, OP>>;

public:
  Solver() = default;
  ~Solver() = default;
  
private:

  template<svdType envalue,
	   typename std::enable_if<
	     envalue==svdType::truncated
	     >::type * = nullptr>
  void computeImpl(matrix_type & A, sc_t tol, int t){
    tol_ = tol;
    //    auto m = A.globalRows();
    auto n = A.globalCols();
    // truncated svd cannot have more than # of cols
    assert(t <= n);

    nU_ = t;
    nV_ = t;
    
    //// comp m left-sing vec as solution of A*A^T U = lambda U
    computeSingVectors(A, t, 'L');
    storeSingVectors('L');

    //// comp n right-sing vec as solution of A^T*A U = lambda U
    computeSingVectors(A, t, 'R');
    storeSingVectors('R');

    computeSingValues(t);
    storeSingValues();
    
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
  
  const sval_t & singularValuesImpl() const{
    return *sval_;
  }//end method
  //-----------------------------------------------------

  
private:
  
  void computeSingVectors(matrix_type & A, int N, char which){
    assert(A.isFillingCompleted());

    //if not already defined, create the target operator for eigensolve
    // build A*A^T or A^T*A
    auto Aout = buildMatForEigSolve(A, which);
    auto Afin = Teuchos::rcpFromRef(*Aout.data() );

    // solve eigenproblem
    solveEigProblem(Afin, N);
    
  }//end method
  //-----------------------------------------------------
  
  void storeSingVectors(char which)
  {
    const Anasazi::Eigensolution<sc_t,MV> & sol = eigP_->getSolution();
    Teuchos::RCP<MV> EV = sol.Evecs;
    if (which == 'L') lsv_ = std::make_shared<lsv_t>(*EV);
    if (which == 'R') rsv_ = std::make_shared<rsv_t>(*EV);

  }//end method
  //-----------------------------------------------------

  auto buildMatForEigSolve(matrix_type & A, char which)
  {
    if (which == 'L'){
      //auto AAT = core::matrixMatrixProduct(A, ASWT, false, false, true);
      auto AAT = core::matrixMatrixProduct(A, A, false, true, true); 
      assert( AAT.isFillingCompleted() );
      return AAT;
   }
   else{
     // auto ASWT = core::transpose(A);
     // auto ATA = core::matrixMatrixProduct(ASWT, A, false, false, true);
     auto ATA = core::matrixMatrixProduct(A, A, true, false, true);
     assert( ATA.isFillingCompleted() );
     return ATA;
    }
    
  }//end method
  //-----------------------------------------------------

      
  void computeSingValues(int N)
  {
    const Anasazi::Eigensolution<sc_t,MV> & sol = eigP_->getSolution();
    std::vector<Anasazi::Value<sc_t> > evals = sol.Evals;
    assert((int)evals.size() == N);
    // compute singular values as sqrt of eigenvalues
    std::transform(evals.cbegin(), evals.cend(),
		   std::back_inserter(eigVals_),
		   [](const Anasazi::Value<sc_t> & value) {
		     return std::sqrt(value.realpart);
		   });
    // //print
    // for (auto it : eigVals_)
    // 	std::cout << it << " "  << std::endl;
  }//end method
  //-----------------------------------------------------

  
  void solveEigProblem(Teuchos::RCP<Epetra_CrsMatrix> & Arcp, int nev)
  {
    //Arcp->Print(std::cout);
    //get domain map of operator
    auto & dommap = Arcp->OperatorDomainMap();
    //dommap.Print(std::cout);

    // create init vectors
    Teuchos::RCP<MV> ivec = Teuchos::rcp(new MV(dommap, blockSize));
    ivec->Random();
    
    // // Construct the Preconditioner
    // Teuchos::RCP<Ifpack_Preconditioner> prec;
    // Teuchos::RCP<Epetra_Operator> PrecOp;
    // Ifpack precFactory;
    // // Set up Ifpack to use incomplete Cholesky with thresholding on
    // // each MPI process and additive Schwarz domain decomposition
    // // across MPI processes.  See Ifpack's documentation for details.
    // std::string precType = "IC stand-alone";
    // int overlapLevel = 0;
    // prec = Teuchos::rcp(precFactory.Create(precType, Arcp.get(), overlapLevel));
    // // parameters for preconditioner
    // Teuchos::ParameterList precParams;
    // precParams.set("fact: drop tolerance",static_cast<sc_t>(1e-4));
    // precParams.set("fact: level-of-fill",0);
    // prec->SetParameters(precParams);
    // prec->Initialize();
    // prec->Compute();
    // //printer.stream(Errors) << " done." << endl;
    // // encapsulate this preconditioner into a IFPACKPrecOp class
    // PrecOp = Teuchos::rcp (new Epetra_InvOperator (&*prec));
    
    // Create the eigenproblem. 
    //Teuchos::RCP<Anasazi::BasicEigenproblem<sc_t, MV, OP>>
    eigP_ = Teuchos::rcp(new Anasazi::BasicEigenproblem<sc_t,MV,OP>(Arcp,ivec));
    eigP_->setNEV(nev);
    eigP_->setHermitian(true);
    //    eigP_->setPrec(PrecOp);
    ///finishing passing it information.
    assert(eigP_->setProblem());
  
    // Create the solver manager and solve
    Teuchos::ParameterList PL;
    PL.set("Which", "LM"); // ordering of eigenvalues, lowest to largest

    PL.set("Block Size", 1);
    PL.set( "Num Blocks", 12);
    PL.set( "Maximum Restarts", 100 );
    PL.set ("Convergence Tolerance", tol_);
    PL.set ("Maximum Iterations", maxIters);
    // PL.set( "Sort Manager", MySort );
    //    PL.set("Dynamic Extra NEV",false);
    //PL.set( "Full Ortho", true );
    //PL.set( "Use Locking", true );
    const int verbosity = Anasazi::Errors +
      Anasazi::FinalSummary +
      Anasazi::Warnings;
    PL.set ("Verbosity", verbosity);

    // solver object
    Anasazi::BlockKrylovSchurSolMgr<double, MV, OP> solverMan(eigP_, PL);
    //Anasazi::BlockDavidsonSolMgr<sc_t, MV, OP> solverMan(eigP_, PL);
    //    Anasazi::LOBPCGSolMgr<sc_t, MV, OP> solverMan(eigP_, PL);

    auto returnCode = solverMan.solve();
    assert(returnCode == Anasazi::Converged );
    
  }// end method
  //-----------------------------------

  
    
  template <typename T = sval_t,
	    typename std::enable_if<
	      std::is_same<T,
		core::Matrix<Epetra_CrsMatrix>
		>::value
	      >::type * = nullptr>
  void storeSingValues()
  {
    auto & comm = lsv_->commCRef();

    Epetra_Map rmap(nU_, 0, comm);
    Epetra_Map cmap(nV_, 0, comm);

    using GO_t = typename core::details::traits<sval_t>::global_ordinal_t;
    // create matrix to store diagonal
    sval_ = std::make_shared<sval_t>(rmap, 1);
    // fill with 1.
    auto myN = rmap.NumMyElements();
    std::vector<GO_t> mygid(myN);
    rmap.MyGlobalElements( mygid.data() );

    std::vector<sc_t> vals(1);
    std::vector<GO_t> inds(1);
    for (auto it : mygid){
      vals[0] = 0.0;
      if (it < static_cast<GO_t>(eigVals_.size()) )
	vals[0] = eigVals_[it];
	
      inds[0] = it;
      sval_->insertGlobalValues(it, 1, vals.data(), inds.data());
    }
    sval_->fillingIsCompleted(cmap, rmap);
    
  }//end method
  //-----------------------------------------------------


  template <typename T = sval_t,
	    typename std::enable_if<
	      std::is_same<T,
		core::Vector<Epetra_Vector>
		>::value
	      >::type * = nullptr>
  void storeSingValues()
  {
    auto & comm = lsv_->commCRef();
    using GO_t = typename core::details::traits<sval_t>::global_ordinal_t;

    // map 
    Epetra_Map mymap( std::min(nU_, nV_), 0, comm);
    // vector object
    sval_ = std::make_shared<sval_t>(mymap);
    // fill
    auto myN = mymap.NumMyElements();
    std::vector<GO_t> mygid(myN);
    mymap.MyGlobalElements( mygid.data() );
    std::vector<sc_t> vals(1);
    std::vector<GO_t> inds(1);
    for (auto it : mygid){
      vals[0] = eigVals_[it];
      inds[0] = it;
      sval_->replaceGlobalValues(1, inds.data(), vals.data());
    }
    
  }//end method
  //-----------------------------------------------------

   
private:
  friend SolverBase< Solver<matrix_type, lsv_type,
			    rsv_type, sval_type> >;

private:
  sc_t tol_;
  const int maxIters = 500;

  // not sure what this is... i cannot find any documentation
  const int blockSize = 5;
  
  std::shared_ptr<lsv_t> lsv_;
  std::shared_ptr<rsv_t> rsv_;

  std::shared_ptr<sval_t> sval_;
  std::vector<sc_t> eigVals_;

  eigprob_t eigP_;
  int nU_;
  int nV_;

};//end class
  
}//end namespace svd
#endif 
