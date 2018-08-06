
#include "CORE_ALL"
#include "ODE_ALL"
#include "SOLVERS_EXP"
#include "apps_burgers1d_epetra.hpp"
//#include "observer.hpp"
#include "experimental/rom_galerkin_explicit_residual_policy.hpp"
// #include "experimental/rom_galerkin_implicit_jacobian_policy.hpp"
// #include "experimental/basis_operator/rom_basis_operator_default.hpp"
// #include "experimental/sampling_operator/rom_sampling_operator_identity.hpp"
// #include "experimental/weighting_operator/rom_weighting_operator_identity.hpp"

struct mysizer{
 using state_t = core::Vector<apps::Burgers1dEpetra::state_type>;
 static size_t getSize(state_t & obj){
   return obj.globalSize();
 };
  static void matchSize(const state_t & src, state_t & obj){
    //obj.resize(src.size());
 };
};

// template<typename T>
// void printSol(std::string mess, const T & y){
//   std::cout << mess << std::endl;
//   for (int i=0; i<y.size(); ++i)
//     std::cout << std::setprecision(14) << y[i]  << " ";
//   std::cout << std::endl;
// }


template <typename matrix_type,
	  typename vector_type,
	  typename std::enable_if<
	    core::details::traits<matrix_type>::isMatrix==1 &&
	    core::details::traits<matrix_type>::isEpetra==1 &&
	    core::details::traits<matrix_type>::isDense==1 &&
	    core::details::traits<vector_type>::isVector==1 &&
	    core::details::traits<vector_type>::isEpetra==1
	    >::type * = nullptr>
auto matrixVectorProduct(const matrix_type & A,
			 const vector_type & b)
{
  assert( A.globalCols() == b.globalSize() );
  vector_type C( A.getDataMap() );
  C.setZero();
  return C;  
}
//------------------------------------------------------

template <typename mat_type,
	  typename std::enable_if<
	    core::details::traits<mat_type>::isEpetra &&
	    core::details::traits<mat_type>::isDense
	    >::type * = nullptr
	  >
auto pseudoInverse(const mat_type & A)
{
  const int nR = A.globalRows();
  const int nC = A.globalCols();
  Epetra_Map newMap(nC, 0, A.commCRef());
  mat_type C(newMap, nR);
  C.setZero();
  return C;
}
//------------------------------------------------------

// this is the operator: (P phi)^+ P 
template <typename phi_t, typename P_t>
class WeightOperator{
public:
  WeightOperator(phi_t const & phi, P_t const & pop)
    : phi_(&phi), pop_(&pop)
  {
    //fullOp_ = (pop*phi_)^+ * pop;
    // auto a1 = phi_->rightMultiply(pop);
    // auto a2 = a1.pseudoInverse();
    // auto a3 = pop.rightMultiply(a2);
  }
  ~WeightOperator() = default;

  template <typename vecT,
  	    typename std::enable_if<
  	      core::details::traits<vecT>::isVector
  	      >::type * = nullptr
  	    >
  auto apply(const vecT & b)
  {
    // P phi
    auto A1 = core::matrixMatrixProduct(*pop_, *phi_);
    std::cout << "P phi: "
	      << A1.globalRows() << " "
	      << A1.globalCols() << std::endl;

    // (P phi)^+ => dense matrix
    auto A2 = pseudoInverse(A1);
    std::cout << "(P phi)^+ :"
	      << A2.globalRows() << " "
	      << A2.globalCols() << std::endl;

    // P b => vector
    auto b1 = core::matrixVectorProduct(*pop_, b);
    std::cout << "P b :"
	      << b1.globalSize() << std::endl;

    // (P phi)^+ P b => vector
    return matrixVectorProduct(A2, b1);
  }
  
private:
  phi_t const * phi_;
  P_t const * pop_;
};
//-------------------------------------------

template <typename T>
void fillColloc(T & A, Epetra_Map map){
  int myNR = map.NumMyElements();
  std::vector<int> mygid(myNR);
  map.MyGlobalElements( mygid.data() );
  std::array<double,1> vals({1.});
  std::array<int,1> colind;
  for (auto const & it : mygid){
    if(it == 0){
      colind = {1};
      A.insertGlobalValues(it, 1, vals.data(), colind.data());
    }
    if(it == 3){
      colind = {5};
      A.insertGlobalValues(it, 1, vals.data(), colind.data());
    }
  }
}


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
  using space_res_t = core::Vector<native_space_residual_type>;
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
  auto & r0n = appObj.getInitialResidual();
  auto & yFOMMap = appObj.getDataMap();
    
  //-------------------------------
  // wrap with core structures
  state_t y0FOM(y0n);
  space_res_t r0FOM(r0n);
    
  //---------------------------------
  // for anasazi, the eigenvectors are stored in a Multivector
  // but for us, they need to be see as a 2d matrix so we wrap it
  // with a core::Matrix<>. This is important. Maybe we
  // should make sure at compile time the phi operator is a core::matrix.
  // for now, the basis vectors have the same row map as states I assume.
  const int numBasisVecs = 12;
  using phi_type = core::Matrix<Epetra_MultiVector>;
  phi_type phi(yFOMMap, numBasisVecs);
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
  fillColloc(P, Pmap);
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
  using A_type = WeightOperator<phi_type, P_type>;
  A_type A(phi, P);

  //---------------------------------------
  // project initial condition phi^T * y0
  auto y = matrixVectorProduct(phiT, y0FOM);
  std::cout << y.globalSize() << std::endl;
  auto res = A.apply(r0FOM);
  std::cout << "res :"
	    << res.globalSize() << std::endl;
  
  // // // residual and jacob policies
  // using res_pol_t = rom::exp::romGalerkinExplicitResidualPolicy<
  //   state_t, space_res_t, model_eval_t, mysizer, phi_type, A_type>;
  // res_pol_t resObj(y, r0FOM, phiT, A);

  // // stepper
  // using stepper_t = ode::ExplicitEulerStepper<
  //   state_t, space_res_t, scalar_t, target_app_t, mysizer>;
  // stepper_t stepperObj(appObj, y, ...);

  // // integrate in time 
  // snapshot_collector collObj;
  // scalar_t final_t = 35;
  // scalar_t dt = 0.01;
  // ode::integrateNSteps(stepperObj, y, 0.0, dt, 3500, collObj);
  
  
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
