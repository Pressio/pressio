#include <iomanip>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <type_traits>
#include <cmath>
#include <fstream>
#include <cassert>
#include <memory>
#include <map>
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_IntVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Time.h"
#include "AztecOO.h"
#include "ml_include.h"
#include "Epetra_LinearProblem.h"
#include "ml_MultiLevelOperator.h"
#include "ml_epetra_utils.h"


const double PI = 3.14159265358979323846;
enum neighTag {left,right,top,bottom};
const double Re_ = 1000.;

using veci = std::vector<int>;
using vecd = std::vector<double>;

template<typename T>
using rcp = std::shared_ptr<T>;

// N = tot number , k = enumerates proc, 0,1,2,3,...
int proc_kToi(int N, int k){ return k/N; }
int proc_kToj(int N, int k){ return k % N; }
void proc_kToij(int & i, int & j, int N, int k){
  j = proc_kToj(N,k);
  i = proc_kToi(N,k);
}
int proc_ijkTok(int iin, int jin, int N){
  return iin*N + jin;
}


int grid_locijToGlobk(int procRank, int li, int lj, int lN){
  return procRank*lN*lN + li*lN + lj;
}
int grid_globkToglobi(int lN, int procRank,
		      int procRanki, int procRankj, int k){
  int localLinIndex = k-procRank*lN*lN;
  return localLinIndex/lN + procRanki*lN;
}
int grid_globkToglobj(int lN, int procRank,
		      int procRanki, int procRankj, int k){
  int localLinIndex = k-procRank*lN*lN;
  return localLinIndex % lN + procRankj*lN;
}

double vorticity(double x, double y, double t){
  double res = -8.0 * PI*PI * std::sin(2.*PI*x)*std::sin(2.*PI*y);
  // std::cout << x << " " << y << " " << std::sin(2.*PI*x) << " "
  // 	    << std::sin(2.*PI*y) << " " << res << std::endl;
  return res;
  //return 2.0*std::sin(q*x)*std::sin(q*y)*std::exp(-2.*q*q*t/Re_);
}
double stream(double x, double y, double t){
  return std::sin(2.*PI*x)*std::sin(2.*PI*y);
  //  return std::sin(q*x)*std::sin(q*y)*std::exp(-2.*q*q*t/Re_);
}



void fillLapla(int lNsq_, rcp<Epetra_CrsMatrix> lapMat_,
	       rcp<Epetra_Map> map_,
	       int gN_, int lN_, int myID_,
	       int myIDi_, int myIDj_,
	       int nPi, int nPj, const std::map<neighTag,int> & neighbors,
	       double dx_)
{  
  veci MyGlobalElements( lNsq_ );
  map_->MyGlobalElements( MyGlobalElements.data() );

  veci Indices(4);
  vecd ValueDiag(1);
  veci IndicesDiag(1);
  double diagVal = 4.0;
  vecd Values{-1, -1, -1, -1};

  int westID, eastID, northID, southID, diagID;
  int localCount = 0;
  for (int li=0; li<lN_; li++){
    for (int lj=0; lj<lN_; lj++)
      {
	int thisGID = MyGlobalElements[localCount]; 
	int gi = grid_globkToglobi(lN_, myID_, myIDi_, myIDj_, thisGID);
	int gj = grid_globkToglobj(lN_, myID_, myIDi_, myIDj_, thisGID);


	if (thisGID != 0 ){
	  
	  westID = (gj==0) ?
	    grid_locijToGlobk(neighbors.at(neighTag::left), li, lN_-1, lN_)
	    : (lj==0) ? thisGID-lNsq_+lN_-1 : thisGID-1;
	  eastID = (gj==gN_-1) ?
	    grid_locijToGlobk(neighbors.at(neighTag::right), li, 0, lN_)
	    : (lj==lN_-1) ? thisGID+lNsq_-lN_+1 : thisGID+1;

	  northID = (gi==0) ?
	    grid_locijToGlobk(neighbors.at(neighTag::top), lN_-1, lj, lN_) :
	    (li==0) ?
	    grid_locijToGlobk(neighbors.at(neighTag::top), lN_-1, lj, lN_) : thisGID-lN_;
	  
	  southID = (gi==gN_-1) ? grid_locijToGlobk(neighbors.at(neighTag::bottom), 0, lj, lN_) :
	    (li==lN_-1) ?
	    grid_locijToGlobk(neighbors.at(neighTag::bottom), 0, lj, lN_) : thisGID+lN_;

	  if ((myID_==1 && li==lN_-1) || (myID_==6 && li==lN_-1))
	    std::cout << " gi = " << gi
		      << " gj = " << gj
		      << " gID = " << thisGID
		      << " west = " << westID << " " 
		      << " east = " << eastID << " " 
		      << " north = " << northID << " " 
		      << " south = " << southID << " " 
		      << std::endl;
	  
	  Indices[0] = westID; 
	  Indices[1] = eastID;
	  Indices[2] = northID;
	  Indices[3] = southID;
	  //off diagonal
	  lapMat_->InsertGlobalValues(thisGID, 4, Values.data(), Indices.data() );
	  // add diagonal
	  ValueDiag[0] = diagVal;
	  IndicesDiag[0] = thisGID;
	  lapMat_->InsertGlobalValues( thisGID, 1, ValueDiag.data(), IndicesDiag.data() );
	}
	else
	{
	  ValueDiag[0] = 1.0;
	  IndicesDiag[0] = thisGID;
	  lapMat_->InsertGlobalValues( thisGID, 1, ValueDiag.data(), IndicesDiag.data() );
	}
	
	localCount++;
      }//end j
  }//end i

  lapMat_->FillComplete();
  //  lapMat_->Print(std::cout);
}//edn method 



void solvePoisson(rcp<Epetra_Vector> vort_,
		  rcp<Epetra_Vector> stream_,
		  double dx_,
		  int lNsq,
		  rcp<Epetra_CrsMatrix> lapMat_,
		  Epetra_MpiComm & ecomm_,
		  rcp<Epetra_Map> myMap_,
		  rcp<Epetra_Vector> xx_,
		  rcp<Epetra_Vector> yy_,
		  int myID_)
{
  for (int i=0; i<lNsq; i++){
    (*vort_)[i] *= -(dx_*dx_);
    (*stream_)[i] = 0.0; //stream((*xx_)[i],(*yy_)[i],0.0);
  }      
  if (myID_==0) {
    (*vort_)[0] = stream((*xx_)[0],(*yy_)[0],0.0);
  }
  
  Epetra_LinearProblem Problem( lapMat_.get(), stream_.get(), vort_.get() );
  AztecOO Solver(Problem);
  Solver.SetAztecOption( AZ_solver, AZ_GMRESR);
  Solver.SetAztecOption( AZ_conv, AZ_r0);
  //Solver.SetAztecOption( AZ_conv, AZ_noscaled);
  // Solver.SetAztecOption( AZ_output, AZ_all );
  // Solver.SetAztecOption( AZ_precond, AZ_Jacobi );
  // Solver.SetAztecOption( AZ_poly_ord, 10 );

  // ================= MultiLevelOperator SECTION ========================
  int nLevels = 10;            // maximum number of levels
  int maxMgLevels = 6;         //
  ML_Set_PrintLevel(10);       // print level (0 silent, 10 verbose)
  ML* ml_handle;               // container of all ML' data
  ML_Create(&ml_handle, maxMgLevels);
  // convert to epetra matrix, put finest matrix into
  // position maxMgLevels - 1 of the hierarchy. NOTE: the matrix
  // is only wrapped (that is, a suitable getrow() function is used),
  // so data in the linear system matrix are NOT replicated.
  EpetraMatrix2MLMatrix(ml_handle, maxMgLevels-1, lapMat_.get() );
  // create an Aggregate object; this will contain information
  // about the aggregation process for each level
  ML_Aggregate *agg_object;
  ML_Aggregate_Create(&agg_object);
  // select coarsening scheme.
  ML_Aggregate_Set_CoarsenScheme_Uncoupled(agg_object);
  // generate the hierarchy. We decided to use decreasing ordering;
  // one can also use ML_INCREASING (in this case, you need to replace
  // maxMgLevels-1 with 0 in call to EpetraMatrix2MLMatrix())
  nLevels = ML_Gen_MGHierarchy_UsingAggregation(ml_handle, maxMgLevels-1,
  						ML_DECREASING, agg_object);
  // // define the ID of the coarsest level
  int coarsestLevel = maxMgLevels - nLevels;
  // // set up some smoothers. 
  int nits = 1;
  for (int level = maxMgLevels-1; level > coarsestLevel; level--){
    ML_Gen_Smoother_Cheby(ml_handle, level, ML_BOTH, 30., 3);
  } 
  // simple coarse solver. You may want to use Amesos to access
  // to a large variety of direct solvers, serial and parallel
  ML_Gen_Smoother_GaussSeidel(ml_handle, coarsestLevel, ML_BOTH,
  			      nits, ML_DEFAULT);
  // generate the solver
  ML_Gen_Solver(ml_handle, ML_MGV, maxMgLevels-1, coarsestLevel);
  // wrap ML_Operator into Epetra_Operator
  ML_Epetra::MultiLevelOperator  MLPrec(ml_handle, ecomm_, *myMap_, *myMap_);
  // ========== End of MultiLevelOperator SECTION ========================

  // //set this operator as preconditioner for AztecOO
  Solver.SetPrecOperator(&MLPrec);

  // solve
  Solver.Iterate(5000, 1e-13);
  ML_Aggregate_Destroy(&agg_object);
  ML_Destroy(&ml_handle);
  
  //Solver.Iterate(5000, 1e-6);
  // vecd status( AZ_STATUS_SIZE );
  // Solver.GetAllAztecStatus( status.data() );
  // double AztecOOExitStatus = status.operator[](1);
  // std::cout << "AztecOOExitStatus = " << AztecOOExitStatus
  // 	      << " " << AZ_normal << std::endl;
  
}





int main(int argc, char *argv[])
{
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);  

  // # of cells along each axis
  int gN_ = 4;

  // proc info
  /*
	   0    myIDj_     nPj_
	   -------------------
       0   |       |         |
	   |   0   |    1    |
	   |       |         |
      i_   -------------------
           |       |         |
	   |   2   |    3    |
      nPi_ |       |         |
	   -------------------
  */
  
  int myID_ = Comm.MyPID();
  int myIDi_, myIDj_;
  int nProc_ = Comm.NumProc();
  int nPi_ = std::sqrt<int>(nProc_);
  int nPj_ = nPi_;
  proc_kToij(myIDi_, myIDj_, nPj_, myID_);
  
  // grid info
  const double L_ = 1.0;
  double dx_ = L_/static_cast<int>(gN_);
  int lN_ = static_cast<int>(gN_)/static_cast<int>(nPi_);
  int lNsq_ = lN_*lN_;
  int totN_ = gN_*gN_;

  // neighbors
  std::map<neighTag,int> neighbors;
  int leftRank = myIDj_==0 ? proc_ijkTok(myIDi_, nPj_-1, nPj_) : myID_-1;
  int rightRank = myIDj_==nPj_-1 ? proc_ijkTok(myIDi_, 0, nPj_) : myID_+1;
  int topRank = myIDi_==0 ? proc_ijkTok(nPi_-1, myIDj_, nPj_) : myID_-nPj_;
  int bottomRank = myIDi_==nPi_-1 ? proc_ijkTok(0, myIDj_, nPj_) : myID_+nPj_;
  neighbors.insert( {neighTag::left,  leftRank} );
  neighbors.insert( {neighTag::right, rightRank} );
  neighbors.insert( {neighTag::top,   topRank} );
  neighbors.insert( {neighTag::bottom,bottomRank} );
  
  std::cout << " r = " << myID_
	    << " nPi_ = " << nPi_ 
	    << " nPj_ = " << nPj_ 
	    << " myIDi_= " << myIDi_
	    << " myIDj_= " << myIDj_
	    << " l,r,t,b = "
	    << leftRank << " " << rightRank << " "
	    << topRank << " " << bottomRank
	    << " lN_ = " << lN_
	    << " dx_ = " << dx_
	    << std::endl;
    
  //---------------------------------
  // map 
  //---------------------------------
  rcp<Epetra_Map> myMap_;
  {
    int NumMyElements;
    veci MyGlobElements;      
    NumMyElements = 0;
    for(int i1=0; i1<lN_; ++i1){
      for(int j1=0; j1<lN_; ++j1){
  	int gID = grid_locijToGlobk(myID_, i1, j1, lN_);
  	MyGlobElements.emplace_back(gID);
  	NumMyElements++;
      }
    }
    myMap_ = std::make_shared<Epetra_Map>(-1, (int) NumMyElements,
  					  MyGlobElements.data(),
  					  0, Comm ) ;
    myMap_->Print(std::cout);
  }//end create map
  
  //---------------------------------
  // fields
  //---------------------------------
  rcp<Epetra_Vector> vort_;
  rcp<Epetra_Vector> stream_;
  rcp<Epetra_Vector> xx_;
  rcp<Epetra_Vector> yy_;
  
  xx_ = std::make_shared<Epetra_Vector>(*myMap_);
  yy_ = std::make_shared<Epetra_Vector>(*myMap_);
  int NumMyElem = myMap_->NumMyElements();
  veci myGel(NumMyElem);
  myMap_->MyGlobalElements(myGel.data());
  int k = 0;
  for (auto const & it : myGel){
    double x = dx_*0.5 + dx_ * grid_globkToglobj(lN_, myID_, myIDi_, myIDj_, it);
    double y = L_ - dx_*0.5 - dx_ * grid_globkToglobi(lN_, myID_, myIDi_, myIDj_, it);
    (*xx_)[k] = x;
    (*yy_)[k] = y;
    k++;
  }    

  // initialize vorticity
  vort_ = std::make_shared<Epetra_Vector>(*myMap_);
  stream_ = std::make_shared<Epetra_Vector>(*myMap_);
  for (int i=0; i<lNsq_; i++){
    (*vort_)[i] = vorticity( (*xx_)[i], (*yy_)[i], 0.0 );
    (*stream_)[i] = 0.0;
    std::cout << i << " "
    	      <<  (*xx_)[i] << " "
    	      <<  (*yy_)[i] << " "
    	      << (*vort_)[i] << " "
    	      << (*stream_)[i] << " "      
    	      << std::endl;
  }

  // poisson matrix
  rcp<Epetra_CrsMatrix> lapMat_;
  int numOfNonZerosEntries = 5;
  lapMat_ = std::make_shared<Epetra_CrsMatrix>(Copy, *myMap_, numOfNonZerosEntries );
  fillLapla(lNsq_, lapMat_, myMap_, gN_, lN_, myID_,
  	    myIDi_, myIDj_, nPi_, nPj_, neighbors, dx_);

  solvePoisson(vort_, stream_, dx_, lNsq_, lapMat_, Comm, myMap_, xx_, yy_, myID_);
  //  stream_->Print(std::cout);
  if (myID_==0)
    {
      double err = 0.0;
      for (int i=0; i<lNsq_; i++){
  	double newv = std::abs( (*stream_)[i] - stream((*xx_)[i],(*yy_)[i],0.0) );
  	err += newv*newv;
  	std::cout << std::setprecision(14)
  		  << (*xx_)[i] << " "
  		  << (*yy_)[i] << " "
  		  << (*stream_)[i] << " "
  		  << stream((*xx_)[i],(*yy_)[i],0.0)
  		  << std::endl;    
      }
      std::cout << "err = " << std::sqrt(err/(gN_*gN_)) << std::endl;
    }
    
  MPI_Finalize() ;  
  return 0;
}



// template <typename state_type, typename oapp>
// class GP()
// {
// private:
//   state_type myY_;
//   oapp * appPtr_;
    
// public:
//   GP(state_type & src, ...)
//     : myY(src), ... {}
//   ~GP(){}
  
//   void operator()()
//   {
//     (*oappPtr)()( V' * y )
//   }

//   void run()
//   {
//     ode::eulerStepper<state_t,state_t,double,stateResizer> myStepper;
//     ode::integrateNSteps(myStepper, *this, myY, snColl, 0.0, 0.0035, 10000);    
//   }
// };




  // bool success = true;
  // std::cout << std::boolalpha << success;

  // MPI_Init(&argc,&argv);
  // int rank; // My process ID
  // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // Epetra_MpiComm Comm(MPI_COMM_WORLD);

  // int MyPID = Comm.MyPID();
  // int NumProc = Comm.NumProc();

  // int NumMyElements = 10000;
  // int NumMyElements1 = NumMyElements; // Needed for localmap
  // int NumGlobalElements = NumMyElements*NumProc+EPETRA_MIN(NumProc,3);
  // if (MyPID < 3) NumMyElements++;
  // int IndexBase = 0;
  // int ElementSize = 7;

  // Epetra_LocalMap *LocalMap = new Epetra_LocalMap(NumMyElements1, IndexBase,Comm);
  // Epetra_Vector A(*LocalMap);

  // // epetramock::evector obj;
  // //obj.print();
 
  //   MPI_Finalize();
