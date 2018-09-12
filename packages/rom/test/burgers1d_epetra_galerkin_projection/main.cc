
#include "CORE_ALL"
#include "ODE_ALL"
#include "SOLVERS_ALL"
#include "apps_burgers1d_epetra.hpp"
#include "experimental/rom_galerkin_explicit_residual_policy.hpp"
#include "experimental/rom_multi_vector_operator.hpp"
#include "observer.hpp"


int main(int argc, char *argv[])
{
  //-------------------------------
  // define native app types
  using app_t = apps::Burgers1dEpetra;
  using app_state_t = app_t::state_type;
  using app_space_res_t = app_t::space_residual_type;
  using scalar_t = app_t::scalar_type;

  //-------------------------------
  // MPI init
  MPI_Init(&argc,&argv);
  int rank; // My process ID
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);

  //-------------------------------
  // create app object
  int numCell = 20; // # fv cells
  std::vector<double> mu({5.0, 0.02, 0.02});
  app_t appObj(mu, numCell, &Comm);
  appObj.setup();
  auto & y0n = appObj.getInitialState();
  auto & r0n = appObj.getInitialResidual();  
  /// wrap with core structures
  using app_state_w_t = core::Vector<app_state_t>;
  using app_space_res_w_t = core::Vector<app_space_res_t>;
  app_state_w_t y0FOM(y0n);
  app_space_res_w_t r0FOM(r0n);

  //-------------------------------
  // types for ode
  using ode_state_t = core::Vector<app_state_t>;
  using ode_res_t = ode_state_t;
  ode_state_t y(y0n);
  ode_res_t r(r0n);
  // stepper wants the types for doing time integration,
  // which might be different than the native ones
  using stepper_t = ode::ExplicitEulerStepper<
    ode_state_t, ode_res_t, app_t>;
  stepper_t stepperObj(appObj, y, r);
    
  //----------------------------------------
  // integrate in time and collect snapshots
  scalar_t fint = 35;
  scalar_t dt = 0.01;
  auto totN = static_cast<unsigned int>(fint/dt);
  snaps_collector obsO(10, appObj.getDataMap(), dt);
  ode::integrateNSteps(stepperObj, y, 0.0, dt, totN , obsO);
  //  obsO.printAll();

  //-------------------------------------------------
  // for time being, fake phi, make phi = snapshots
  auto phi = obsO.getSnaps();
  int romSize = phi.globalNumVectors(); // # of snapshots 
  // define the operator wrapping phi
  using phi_op_t = rom::MultiVectorOperator<decltype(phi)>;
  phi_op_t phiOp(phi);
  
  //-------------------------------------------------
  // residual policy needs to know about the APP types
  using res_pol_t = rom::exp::RomGalerkinExplicitResidualPolicy<
    app_state_w_t, app_space_res_w_t, phi_op_t, phi_op_t>;
  res_pol_t resObj(y0FOM, r0FOM, phiOp, phiOp);
  
  /*-------------------------------
    define rom types to use
    this does not have to be the same as natives types, 
    but you can choose what type you want to do ROM */
  using rom_state_t = core::Vector<Eigen::VectorXd>;
  using rom_res_t = rom_state_t;
  rom_state_t yROM(romSize);
  yROM.putScalar(1.0);
  rom_res_t rROM(romSize);

  // stepper needs to know about types for doing time integration,
  // which here are those for ROM
  using rom_stepper_t = ode::ExplicitEulerStepper<
    rom_state_t, rom_res_t, app_t, res_pol_t>;
  rom_stepper_t romStepperObj(appObj, resObj, yROM, rROM);

  /// integrate in time 
  scalar_t dtrom = dt;
  ode::integrateNSteps(romStepperObj, yROM, 0.0, dtrom, 1);

  if (rank==0){
    for (size_t i=0; i<yROM.size(); i++)
      std::cout << yROM[i] << " ";
    std::cout << std::endl;
  }

  MPI_Finalize();
  return 0;
}












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





    
  
  // //-------------------------------
  // // Galerkin ROM
  // //-------------------------------
  // // operators
  // using phiOp_t = rom::experimental::basisOperatorDefault<basis_type>;
  // phiOp_t phiOp(phi);
  // using scale_op_t = rom::experimental::weightingOperatorIdentity<jac_t,phiOp_t>;
  // scale_op_t WOp(phiOp, phi.rows(), phi.cols());

  // // project initial condition
  // size_t redSize = phi.cols();
  // state_t yr(redSize);
  // phiOp.project(y0, yr);
  // // std::cout << "y " << *y.data() << std::endl;
  
  // // residual and jacob policies
  // using res_pol_t = rom::exp::romGalerkinImplicitResidualPolicy<
  //   state_t, residual_t, model_eval_t, scalar_t,
  //   mysizer, phiOp_t, scale_op_t>;
  // res_pol_t resObj(y0, yr, phiOp, WOp);
  // using jac_pol_t = rom::exp::romGalerkinImplicitJacobianPolicy<
  //   state_t, jac_t, model_eval_t, scalar_t, mysizer, phiOp_t, scale_op_t>;
  // jac_pol_t jaObj(y0, yr, phiOp, WOp);

  // //-----------------------------------------------
  // // SOLVERS
  // using lin_solve_t
  //   = solvers::experimental::linearSolver<jac_t, state_t, state_t>;
  // lin_solve_t ls;
  // using nonlin_solve_t
  //   = solvers::experimental::newtonRaphson<state_t, state_t, jac_t, lin_solve_t>;
  // nonlin_solve_t nonls(ls);
  
  // //-----------------------------------------------
  // // stepper
  // using stepper_t = ode::ImplicitEulerStepper<
  //   state_t, residual_t, jac_t, model_eval_t, mysizer,
  //   res_pol_t, jac_pol_t>;
  // stepper_t stepperObj(appObj, resObj, jaObj);

  // //-----------------------------------------------
  // // integrator
  // scalar_t dt = 0.01;
  // scalar_t final_t = dt*1;
  // auto numSteps = static_cast<unsigned int>(final_t/dt);  
  // ode::integrateNSteps(stepperObj, yr, 0.0, dt, numSteps, nonls);

  // //-----------------------------------------------
  // // process final state
  // state_t yrFin(y0);
  // phiOp.leftMultiply(yr, yrFin);
  // // state_t tmp( *phi.data() * (*gg.data()) ) ;
  // printSol("", yrFin);
  // //-------------------------------



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



  // using native_phi_type = Epetra_MultiVector;
  // using my_phi_type = core::MultiVector<native_phi_type>;  


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





// // template<typename T>
// // void printSol(std::string mess, const T & y){
// //   std::cout << mess << std::endl;
// //   for (int i=0; i<y.size(); ++i)
// //     std::cout << std::setprecision(14) << y[i]  << " ";
// //   std::cout << std::endl;
// // }

// template <typename T>
// void fillColloc(T & A, Epetra_Map map){
//   int myNR = map.NumMyElements();
//   std::vector<int> mygid(myNR);
//   map.MyGlobalElements( mygid.data() );
//   std::array<double,1> vals({1.});
//   std::array<int,1> colind;
//   for (auto const & it : mygid){
//     if(it == 0){
//       colind = {1};
//       A.insertGlobalValues(it, 1, vals.data(), colind.data());
//     }
//     if(it == 3){
//       colind = {5};
//       A.insertGlobalValues(it, 1, vals.data(), colind.data());
//     }
//   }
// }
















  // //-------------------------------
  // /* lets say I obtain here the basis somehow
  //    phi type is same as those of the app, so if 
  //    app is trilinos based, phi is also based on trilinos 
  // */
  // auto & yFOMMap = appObj.getDataMap();
  // using phi_type = core::MultiVector<Epetra_MultiVector>;
  // phi_type phi(yFOMMap, 15);
  // using phi_op_t = core::MultiVectorOperator<phi_type>;
  // phi_op_t phiOp(phi);

  // if (rank==0){
  //   //-------------------------------
  //   /* define rom types to use
  //      this does not have to be the same as natives types, 
  //      but you can choose what type you want to do ROM 
  //   */
  //   using rom_state_t = core::Vector<Eigen::VectorXd>;
  //   using rom_res_t = rom_state_t;
  //   using rom_jac_t = core::Matrix<Eigen::MatrixXd>;
  //   rom_state_t yROM;
  //   rom_res_t rROM;
  
  //   // // residual policy needs to know about the app types
  //   using res_pol_t = rom::exp::RomGalerkinExplicitResidualPolicy<
  //     app_state_t, app_space_res_t, phi_op_t, phi_op_t>;
  //   res_pol_t resObj(y0n, r0n, phiOp, phiOp);

  //   // stepper needs to know about types for doing time integration,
  //   // which here are those for ROM
  //   using stepper_t = ode::ExplicitEulerStepper<
  //     rom_state_t, rom_res_t, app_t, res_pol_t>;
  //   stepper_t stepperObj(appObj, resObj, yROM, rROM);

  //   // integrate in time 
  //   scalar_t dt = 0.01;
  //   ode::integrateNSteps(stepperObj, yROM, 0.0, dt, 1);
  // }





  // //---------------------------------
  // // for anasazi, the eigenvectors are stored in a Multivector
  // // but for us, they need to be see as a 2d matrix so we wrap it
  // // with a core::Matrix<>. This is important. Maybe we
  // // should make sure at compile time the phi operator is a core::matrix.
  // // for now, the basis vectors have the same row map as states I assume.
  // const int numBasisVecs = 12;
  // using phi_type = core::Matrix<Epetra_MultiVector>;
  // phi_type phi(yFOMMap, numBasisVecs);
  // auto phiT = core::transpose(phi);
  
  // //---------------------------------
  // // collocation operator
  // const int sampleMeshSize = 18; //has to be > numBasisVecs
  // Epetra_Map Pmap(sampleMeshSize, 0, Comm);
  // using P_type = core::Matrix<Epetra_CrsMatrix>;
  // // 1 nonzero for each row because P contains rows of
  // // the identity matrix of size numCell x nunCell.
  // // So after we select only subset of rows, then P
  // // should have size sampleMesh x numCell
  // P_type P(Pmap, 1);
  // // fill collocation operator
  // fillColloc(P, Pmap);
  // // ...
  // // when we fill complete P, we need to make sure we pass the
  // // domain and range maps so we can do operatiors because P is a rect matrix.
  // // (a) P will leftmultiply phi, so it has to have domain map of phi.
  // // (b) P will leftmultiply the space residual, so it has to have
  // // domain map of the rhs f(...)
  // P.fillingIsCompleted(yFOMMap, Pmap);
  // P.data()->Print(std::cout);
  
  // //---------------------------------
  // // weighting operator
  // using A_type = rom::exp::WeightOperator<phi_type, P_type>;
  // A_type A(phi, P);

  // //---------------------------------------
  // // project initial condition phi^T * y0
  // auto yROM = core::matrixVectorProduct(phiT, y0FOM);
  // std::cout << "yromSz = " << yROM.globalSize() << std::endl;
  // auto rROM = A.apply(r0FOM);
  // std::cout << "rromSz = " << rROM.globalSize() << std::endl;
 
  // // residual policy
  // using res_pol_t = rom::exp::RomGalerkinExplicitResidualPolicy<
  //   state_t, space_res_t, model_eval_t, phi_type, A_type>;
  // res_pol_t resObj(y0FOM, r0FOM, phi, A);

  // // // stepper
  // // using stepper_t = ode::ExplicitEulerStepper<
  // //   state_t, space_res_t, model_eval_t, res_pol_t>;
  // // stepper_t stepperObj(appObj, resObj, yROM, rROM);

  // // // integrate in time 
  // // //snapshot_collector collObj;
  // // //  scalar_t final_t = 35;
  // // scalar_t dt = 0.01;
  // // ode::integrateNSteps(stepperObj, yROM, 0.0, dt, 1);//, collObj);
