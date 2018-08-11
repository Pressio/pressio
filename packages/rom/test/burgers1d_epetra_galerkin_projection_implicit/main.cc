
#include "CORE_ALL"
#include "ODE_ALL"
#include "SOLVERS_EXP"
#include "apps_burgers1d_epetra.hpp"
#include "experimental/rom_galerkin_implicit_residual_policy.hpp"
#include "experimental/rom_galerkin_implicit_jacobian_policy.hpp"
#include "experimental/rom_matrix_pseudo_inverse.hpp"
#include "experimental/rom_operators.hpp"

template <typename T>
void fillColloc(T & A, Epetra_Map map, int sampleMeshSize)
{
  // collocation matrix contains target rows of the identity
  // matrix of those points where to sample
  std::vector<int> targRID {
      0,3,4,5,6,7,18,45,46,47,68,69,81,82,83,84,85,86};
  
  int myNR = map.NumMyElements();
  std::vector<int> mygid(myNR);
  map.MyGlobalElements( mygid.data() );

  std::array<double,1> vals({1.});
  std::array<int,1> colind;
  for (auto const & it : mygid)
  {
    colind = {targRID[it]};
    A.insertGlobalValues(it, 1, vals.data(), colind.data());
  }
};

struct solverFake
{  
  template <typename T, typename T2>
  void solve(T & y, T2 & b){
    
  }
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
  //-------------------------------
  // define native types
  using model_eval_t = apps::Burgers1dEpetra;
  using native_state_t = model_eval_t::state_type;
  using native_space_residual_type = model_eval_t::space_residual_type;
  using native_jac_t = model_eval_t::jacobian_type;
  using scalar_t = model_eval_t::scalar_type;

  //-------------------------------
  // define wrapper types
  using state_t = core::Vector<native_state_t>;
  using res_t = core::Vector<native_space_residual_type>;
  using jac_t = core::Matrix<native_jac_t>;

  //-------------------------------
  // MPI init
  MPI_Init(&argc,&argv);
  int rank; // My process ID
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);

  //-------------------------------
  // create app object
  int numCell = 100; // number of fv cells
  std::vector<double> mu({5.0, 0.02, 0.02});
  model_eval_t appObj(mu, numCell, &Comm);
  appObj.setup();
  auto & y0n = appObj.getInitialState();
  auto & yFOMMap = appObj.getDataMap();
  auto & r0n = appObj.getInitialResidual();
  auto & j0n = appObj.getInitialJacobian();
    
  //-------------------------------
  // wrap with core structures
   state_t y0FOM(y0n);
   res_t r0FOM(r0n);
   jac_t j0FOM(j0n);
    
  //---------------------------------
  // for anasazi, the eigenvectors are stored in a Multivector
  // but for us, they need to be see as a 2d matrix so we wrap it
  // with a core::Matrix<>. This is important. Maybe we
  // should make sure at compile time the phi operator is a core::matrix.
  // for now, the basis vectors have the same row map as states I assume.
  const int numBasisVecs = 12;
  using phi_type = core::Matrix<Epetra_MultiVector>;
  phi_type phi(yFOMMap, numBasisVecs);
  phi.data()->Random();
  phi.data()->Print(std::cout);
  auto phiT = core::transpose(phi);
  
  //---------------------------------
  // collocation operator
  const int sampleMeshSize = 18; //has to be > numBasisVecs
  Epetra_Map Pmap(sampleMeshSize, 0, Comm);
  using P_type = core::Matrix<Epetra_CrsMatrix>;
  // 1 nonzero for each row because P contains rows of
  // the identity matrix of size numCell x nunCell.
  // So after we select only subset of rows, then P
  // should have size sampleMesh x numCell
  P_type P(Pmap, 1);
  // fill collocation operator
  fillColloc(P, Pmap, sampleMeshSize);
  // ...
  // when we fill complete P, we need to make sure we pass the
  // domain and range maps so we can do operatiors because P is a rect matrix.
  // (a) P will leftmultiply phi, so it has to have domain map of phi.
  // (b) P will leftmultiply the space residual, so it has to have
  // domain map of the rhs f(...)
  P.fillingIsCompleted(yFOMMap, Pmap);
  P.data()->Print(std::cout);
  
  //---------------------------------
  // weighting operator
  using A_type = rom::exp::WeightOperator<phi_type, P_type>;
  A_type A(phi, P);

  //---------------------------------------
  // project initial condition phi^T * y0
  //  auto yROM = core::matrixVectorProduct(phiT, y0FOM);
  //  std::cout << "yromSz = " << yROM.globalSize() << std::endl;

  //  j0FOM.data()->Print(std::cout);
  
  // // residual and jacob policies
  // using res_pol_t = rom::exp::RomGalerkinImplicitResidualPolicy<
  //   state_t, res_t, model_eval_t, phi_type, A_type>;
  // res_pol_t resObj(y0FOM, r0FOM, phi, A);

  // using jac_pol_t = rom::exp::RomGalerkinImplicitJacobianPolicy<
  //   state_t, jac_t, model_eval_t, phi_type, A_type>;
  // jac_pol_t jaObj(y0FOM, j0FOM, phi, A);
  
  // // stepper
  // using stepper_t = ode::ImplicitEulerStepper<
  //   state_t, res_t, jac_t, model_eval_t, res_pol_t, jac_pol_t>;
  // stepper_t stepperObj(appObj, resObj, jaObj, yROM);


  // // define solver
  // /* solver needs to know how to construct residual and jacobian objects
  //    so we need to pass them to the solver because solver is the one
  //    holding the residual and jacobian objects.
  //    Need to be careful because solver needs to have the size of the 
  //    REDUCED system, so taking into account the weighting operator.
  //    Also, the type of the residual and jacobain need to be carefully set. 
  //    FOr instance, the type of the jacobian should be inferred from 
  //    the return type of the Apply method of the weighting operator. 
  // */
  // solverFake sobj;
  // // setup

  // //   // integrate in time 
  // // //  scalar_t final_t = 35;
  // // scalar_t dt = 0.01;
  // // ode::integrateNSteps(stepperObj, yROM, 0.0, dt, 1, sobj);
    
  MPI_Finalize();
  return 0;
}






// {
//   // P phi
//   auto Pphi = core::matrixMatrixProduct(P, phi);
//   std::cout << "P phi: "
// 	      << Pphi.globalRows() << " "
// 	      << Pphi.globalCols() << std::endl;
//   //    Pphi.data()->Print(std::cout);

//   // (P phi)^+ => dense matrix
//   auto PphiPI = rom::exp::pseudoInverse(Pphi);
//   PphiPI.data()->Random();
//   std::cout << "(P phi)^+ :"
//   	      << PphiPI.globalRows() << " "
//   	      << PphiPI.globalCols() << std::endl;
//   //    PphiPI.data()->Print(std::cout);

//   // J phi => matrix
//   auto Jphi = core::matrixMatrixProduct(j0FOM, phi);
//   std::cout << "J phi :"
//   	      << Jphi.globalRows() << " "
//   	      << Jphi.globalCols() << std::endl;
    
//   // P Jphi => matrix
//   auto PJphi = core::matrixMatrixProduct(P, Jphi);
//   std::cout << "P Jphi :"
//   	      << PJphi.globalRows() << " "
//   	      << PJphi.globalCols() << std::endl;
    
//   // (P phi)^+ P J phi => matrix
//   auto JJF = core::matrixMatrixProduct(PphiPI, PJphi);
//   std::cout << "JJF :"
//   	      << JJF.globalRows() << " "
//   	      << JJF.globalCols() << std::endl;
//   //    JJF.data()->Print(std::cout);

// }




// Eigen::MatrixXd readPhi(int nr, int nc)
// {
//   Eigen::MatrixXd phi;
//   phi.resize(nr,nc);  
//   std::ifstream source;
//   source.open("bas.txt", std::ios_base::in);
//   std::string line;
//   int row = 0;
//   while (std::getline(source, line) ){
//     //make a stream for the line itself
//     std::istringstream in(line);
//     // tmp variable to store each entry of the file
//     std::vector<std::string> cols(nc);
//     for (int i=0; i<nc; i++){
//       in >> cols[i];
//       phi(row, i) = atof(cols[i].c_str());
//     }
//     row++;
//   }
//   source.close();
//   return phi;
// }//end 

// template <typename container_type>
// class RomOperator{
// public:
//   RomOperator() = default;
//   ~RomOperator() = default;

// private:
//   container_type * A_;
// };
// //-------------------------------------------





  // //-----------------------------------------------
  // // SOLVERS
  // using lin_solve_t
  //   = solvers::experimental::linearSolver<jac_t, state_t, state_t>;
  // lin_solve_t ls;
  // using nonlin_solve_t
  //   = solvers::experimental::newtonRaphson<state_t, state_t, jac_t, lin_solve_t>;
  // nonlin_solve_t nonls(ls);
  


  // using stepper_t = ode::ImplicitEulerStepper<
  //   state_t, residual_t, jac_t, scalar_t, model_eval_t,
  //   scalar_t, mysizer, nonlin_solve_t, res_pol_t, jac_pol_t>;
  // stepper_t stepperObj(appObj, nonls, resObj, jaObj);  
  // snapshot_collector collObj;
  // ode::integrateNSteps(stepperObj, y, 0.0, dt, numSteps, collObj);
  // printSol("", y+y0);
  // //using stepper_t = ode::ExplicitRungeKutta4Stepper<
  // using stepper_t = ode::ExplicitEulerStepper<
  //   state_t, residual_t, scalar_t, model_eval_t,
  //   scalar_t, mysizer>;//, res_pol_t>;
  // stepper_t stepperObj(appObj);//, resObj);
  // // // integration details
  // scalar_t dt = 0.01;
  // scalar_t final_t = 35.;//dt*100;
  // snapshot_collector collObj;
  // ode::integrateNSteps(stepperObj, y, 0.0, dt, final_t/dt, collObj);    
  // printSol("", y);


  // //-------------------------------
  // // collect snapshots usinf FOM
  // //-------------------------------
  // state_t y(y0);
  // snapshot_collector collObj(numCell, numSteps);

  // using stepper_t = ode::ImplicitEulerStepper<
  //   state_t, residual_t, jac_t, scalar_t, model_eval_t,
  //   scalar_t, mysizer, nonlin_solve_t>;
  // stepper_t stepperObj(appObj, nonls);
  
  // ode::integrateNSteps(stepperObj, y, 0.0, dt, numSteps, collObj);
  // printSol("", y);
  // // //  collObj.printAll();

  // //-------------------------------
  // // SVD
  // //-------------------------------
  // Eigen::JacobiSVD<Eigen::MatrixXd> svd(collObj.snapshots_, Eigen::ComputeThinU);
  // auto phi_nat = svd.matrixU();
  // auto phiT_nat = phi_nat;
  // phiT_nat.transposeInPlace();
  // std::cout << "phiSize " << phi_nat.rows() << " " << phi_nat.cols() << std::endl;
  // // std::cout << phi_nat << std::endl;
  // // std::cout  << std::endl;
  // // std::cout << phiT_nat << std::endl;
