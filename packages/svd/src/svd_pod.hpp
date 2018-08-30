
#ifndef ROM_POD_HPP_
#define ROM_POD_HPP_

#include "rom_ConfigDefs.hpp"
#include "CORE_MATRIX"
//#include "AnasaziBlockDavidsonSolMgr.hpp"
//#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include <Teuchos_RCPDecl.hpp>
#include "AnasaziLOBPCGSolMgr.hpp"
#include "AnasaziSimpleLOBPCGSolMgr.hpp"
#include "Epetra_InvOperator.h"

namespace rom{
namespace exp{


template <typename mat_type,
	  typename std::enable_if<
	    core::details::traits<mat_type>::isEpetra &&
	    core::details::traits<mat_type>::isSparse
	    >::type * = nullptr
	  >
auto podImpl(const mat_type & A, int nev,
	     typename core::details::traits<mat_type>::scalar_t tol,
	     int blockSize, int maxIters,
	     bool usePreconditioner = false)
{

  using sc_t = typename core::details::traits<mat_type>::scalar_t;

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

  // //-----------------------------------------------------------------
  // // Construct the Preconditioner
  // RCP<Ifpack_Preconditioner> prec;
  // RCP<Epetra_Operator> PrecOp;
  // if (usePrec) {
  //   printer.stream(Errors) << "Constructing Incomplete Cholesky preconditioner..." << std::flush;
  //   Ifpack precFactory;
  //   // Set up Ifpack to use incomplete Cholesky with thresholding on
  //   // each MPI process and additive Schwarz domain decomposition
  //   // across MPI processes.  See Ifpack's documentation for details.
  //   std::string precType = "IC stand-alone";
  //   int overlapLevel = 0;
  //   prec = rcp (precFactory.Create (precType, K.get (), overlapLevel));
  //   // parameters for preconditioner
  //   Teuchos::ParameterList precParams;
  //   precParams.set("fact: drop tolerance",prec_dropTol);
  //   precParams.set("fact: level-of-fill",prec_lofill);
  //   IFPACK_CHK_ERR(prec->SetParameters(precParams));
  //   IFPACK_CHK_ERR(prec->Initialize());
  //   IFPACK_CHK_ERR(prec->Compute());
  //   //
  //   printer.stream(Errors) << " done." << endl;
  //   // encapsulate this preconditioner into a IFPACKPrecOp class
  //   PrecOp = rcp (new Epetra_InvOperator (&*prec));
  // }
  
  //-----------------------------------------------------------------
  // Create the eigenproblem.  This object holds all the stuff about
  // your problem that Anasazi will see. 
  RCP<Anasazi::BasicEigenproblem<double, MV, OP> > problem =
    rcp (new Anasazi::BasicEigenproblem<double, MV, OP> (Arcp, ivec));
  problem->setNEV (nev);
  problem->setHermitian(true);
  //problem->setPrec (PrecOp);
  // Tell the eigenproblem that you are finishing passing it information.
  assert(problem->setProblem());
  
  //-----------------------------------------------------------------
  // Create the solver manager and solve
  Teuchos::ParameterList PL;
  PL.set ("Which", "LM");
  PL.set ("Block Size", blockSize);
  PL.set ("Maximum Iterations", maxIters);
  PL.set ("Convergence Tolerance", tol);
  const int verbosity = Anasazi::Errors + /*Anasazi::Warnings +*/ Anasazi::FinalSummary;
  PL.set ("Verbosity", verbosity);

  // actual solver object
  Anasazi::SimpleLOBPCGSolMgr<double, MV, OP> solverMan(problem, PL);
  auto returnCode = solverMan.solve();
  assert(returnCode == Anasazi::Converged );

  //-----------------------------------------------------------------
  // get the solution
   const Anasazi::Eigensolution<sc_t,MV> & sol = problem->getSolution();
   const int numCompEV = sol.numVecs;
   std::vector<Anasazi::Value<sc_t> > evals = sol.Evals;
   RCP<MV> evecs = sol.Evecs;
   return core::MultiVector<MV>(*evecs);
   //   evecs->Print(std::cout);
   // for (auto it : evals)
   //   std::cout << it.realpart << " "  << std::endl;  

}
//--------------------------------------------------
  
  
template <typename mat_type,
	  typename std::enable_if<
	    core::details::traits<mat_type>::isEpetra &&
	    core::details::traits<mat_type>::isSparse
	    >::type * = nullptr
	  >
auto pod(mat_type & A, // non-const because otherwise transpose does not work
	 bool isHermitian, 
	 int nev,
	 typename core::details::traits<mat_type>::scalar_t tol = 1e-12,
	 int blockSize = 2, //this is tricky, leave it this way
	 int maxIters = 500,
	 bool usePreconditioner = false)
{

  auto nrA = A.globalRows();
  auto ncA = A.globalCols();
  if ( nrA<=ncA)
    assert( nev <= nrA );
  if ( nrA>ncA)
    assert( nev <= ncA );

  
  if (!A.isFillingCompleted())
	 A.fillingIsCompleted();
	 
  if (isHermitian){
    return podImpl(A, nev, tol, blockSize, maxIters, usePreconditioner);
  }
  else{
    auto ASWT = core::transpose(A);
    auto AAT = core::matrixMatrixProduct(A, ASWT, false, false, true);
    AAT.data()->Print(std::cout);
    return podImpl(AAT, nev, tol, blockSize, maxIters, usePreconditioner);
  }

}
//------------------------------------------------------

  
}//end namespace exp
}//end namespace rom
#endif 













// // Create the Block Davidson eigensolver.
// Anasazi::BlockDavidsonSolMgr<double,MV,OP> anasaziSolver (problem, PL);
// Anasazi::ReturnType returnCode = anasaziSolver.solve ();
// if (returnCode != Anasazi::Converged) {
//   std::cout << "Anasazi eigensolver did not converge." << std::endl;
// }
// // Get the eigenvalues and eigenvectors from the eigenproblem.
// Anasazi::Eigensolution<double,MV> sol = problem->getSolution ();
// // Anasazi returns eigenvalues as Anasazi::Value, so that if
// // Anasazi's Scalar type is real-valued (as it is in this case), but
// // some eigenvalues are complex, you can still access the
// // eigenvalues correctly.  In this case, there are no complex
// // eigenvalues, since the matrix pencil is symmetric.
// std::vector<Anasazi::Value<double> > evals = sol.Evals;
// RCP<MV> evecs = sol.Evecs;
// // Print the results on MPI process 0.
// if (A.commCRef().MyPID() == 0){
//   std::cout << "Solver manager returned "
//        << (returnCode == Anasazi::Converged ? "converged." : "unconverged.")
//        << std::endl << std::endl
//        << "------------------------------------------------------" << std::endl
//        << std::setw(16) << "Eigenvalue"
//        << std::setw(18) << "Direct Residual"
//        << std::endl
//        << "------------------------------------------------------" << std::endl;
//   for (int i=0; i<sol.numVecs; ++i) {
//     std::cout << std::setw(16) << evals[i].realpart
//          << std::endl;
//   }
//   std::cout << "------------------------------------------------------" << std::endl;
// }


// // /////////////////////////////////////////////////////////
// // /////////////////////////////////////////////////////////
// // Create the solver manager
// Anasazi::LOBPCGSolMgr<double, MV, OP> MySolverMan(problem, PL);
// // Solve the problem
// auto returnCode = MySolverMan.solve();

      
// /////////////////////////////////////////////////////////
// /////////////////////////////////////////////////////////
// // Create the Block Krylov-Schur eigensolver.
// Anasazi::BlockKrylovSchurSolMgr<double, MV, OP> MySolverMgr (problem, PL);
// Anasazi::ReturnType returnCode = MySolverMgr.solve ();
// if (returnCode != Anasazi::Converged) {
//   std::cout << "Anasazi eigensolver did not converge." << std::endl;
// }

 
// // Create a set of initial vectors to feed to the eigensolver.
// core::vector<Epetra_MultiVector> ivec(rowmap, blockSize);
// ivec.data()->Random();
    
// const int nR = A.globalRows();
// const int nC = A.globalCols();
// // Epetra_Map newMap(nC, 0, A.commCRef());
// // mat_type C(newMap, nR);
// // C.setZero();
// return C;
// if sz is not None and engy is not None:
//     print('Warning: Specified both sz and engy; sz ignored.')
// if svd is None:
//     from scipy.linalg import svd as svdloc
//     svd = lambda Y: svdloc(Y, full_matrices=False, compute_uv=True)
// U, s, V = svd(X)
// rnk = sum(s>=rnkdef_tol)
// if szbnds is None:
//     szbnds = (1, rnk)
// szbnds = (max(1, szbnds[0]), min(rnk, szbnds[1]))
// if engy is not None:
//     sz = rnk if engy >= 1.0 else np.nonzero(trunc_engy(s)>=engy)[0][0]+1
// if sz is not None:
//     sz = min(max(sz, szbnds[0]), szbnds[1])
// return U[..., :sz], s[:sz], V[..., :sz]

