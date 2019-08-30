
//FRizzi: @Chi, please fix all this 


/* 2 big questions:
	a. what data do you need?
	b. what methods/operations do you need to handle this data?
 */

int main(int argc, char *argv[]){

	//====== step 0: run FOM simulations (classical time marching)? ======//
	/* 
	- run FOM simulations to collect FOM solutions at all time steps which are also snapshots. Then there are 2 ways of creating space-time bases: 

	(a) Classical space-time bases: arrange the solutions at all time steps into one big/long column vector (one long column vector corresponds to one input parameter), then do SVD on these snapshots to create space-time bases, finally re-arrange this long column vector into a matrix of spatial bases at all time instances. This step can be done within pressio or Matlab? This should be done first to verify the solver. 

	(b) Kronecker product space-time bases: tensor product the spatial bases with temporal bases to create space-time bases. This is implemented after (a) above is done. 
	*/

	using fom_t		= pressio::apps::Burgers1dEpetra;
	using scalar_t	= typename fom_t::scalar_type;
	using eig_dyn_vec	= Eigen::Matrix<scalar_t, -1, 1>;
	using lspg_state_t= pressio::containers::Vector<eig_dyn_vec>;

	using decoder_jac_t	= pressio::containers::MultiVector<Epetra_MultiVector>;
  using decoder_t	= pressio::rom::LinearDecoder<decoder_jac_t>;

  std::string checkStr {"PASSED"};

  //-------------------------------
  // MPI init
  MPI_Init(&argc,&argv);
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  if(Comm.NumProc() != 2) return 0;


	//====== step 1: create the app object ======//
	/* we can use any of the time-dependent problems that Francesco has developed. Assume Ns=spatial dof, Nt=temporal dof, nst=number of space-time bases. Set up ROM problem with all these parameters at this step. 

	Maybe include a small function to check the feasibility of the problem based on 3 parameters: Ns, Nt and nst? Namely, we may put some upper bound limits for these 3 parameters.
	
	Data for this step:
	- Nx, Ny, Nt, (Ns=Nx x Ny?), dt, time marching scheme (explicit/implicit, Backward/Forward Euler?) 

	Methods/operations for this step:
	- compute residual, jacobian, source term (i.e., velocity), 
	*/

  //-------------------------------
  // app object
  constexpr int numCell = 20;
  fom_t appObj( {5.0, 0.02, 0.02}, numCell, &Comm);
  appObj.setup();
  auto t0 = static_cast<scalar_t>(0);
  scalar_t dt = 0.01;


	//====== step 2: create the space-time bases ======//
	/*--- step 2.1: create one class to store data: class MultiSpaceTimeBasesData. This class will store data in the MultiVector format (so-called (b) format): 

	[ basis_1__timeStep_1, basis_2__timeStep_1,..., basis_nst__timeStep_1, 
	basis_1__timeStep_2, basis_2__timeStep_2, ..., basis_nst__timeStep_2, 
	... ... 
	basis_1__timeStep_Nt, basis_2__timeStep_Nt, ..., basis_nst__timeStep_Nt ]

	Size of this resulting matrix is (Ns x nst Nt). 
	We will read this data from an input file and may arrange in the above format. Another alternative format will be (so-called (a) format)

	[ basis_1__timeStep_1, basis_1__timeStep_2,..., basis_1__timeStep_Nt, 
	basis_2__timeStep_1, basis_2__timeStep_2, ..., basis_2__timeStep_Nt, 
	... ...
	basis_nst__timeStep_1, basis_nst__timeStep_2, ..., basis_nst__timeStep_Nt ]
	
	Size of this resulting matrix is (Ns x Nt nst). 
	Format (a) here will be considered later when the implementation of format (b) is done.	

	Data for this step: space-time bases collected from step 0 above
	*/

	/*--- step 2.2: one class to store methods: class MultiSpaceTimeBasesMethods. There are 4 methods that we need to define:
	2.2.1 MultiSpaceTimeBasesMethods::MultiplyRomVector
	Calculate U_st = Phi_st * uhat (strictly for Phi_st only)
	Size of Phi_st is Ns x nst Nt, size of uhat is nst, the resulting matrix U_st has size Ns x Nt, i.e., ROM solutions (full size) at all time steps.

	Possible implementation: can implement the product in a way such that we output vector U at time step n_t (i.e., size Ns x 1). 

	Data for this sub-step: uhat, U_st

	2.2.2 MultiSpaceTimeBasesMethods::CalculateJacobianTimesPhiSpaceTime
	Calculate JP = space-time Jac*Phi_st (strictly for Phi_st only)
	Space-time Jacobain comprises spatial Jacobian at all Nt time steps. Size of spatial Jac is Ns x Ns (sparse), size of Phi_st is Ns x nst Nt, size of resulting matrix JP is Ns x nst Nt. 
	Note that JP will have exactly similar size as Phi_st, hence we can use this container to hold JP. (Above we use this container to hold Phi_st.) Look at page 18 in the report to see how we compute each term JP_ij:
	
	JP_ij = J_i_i Phi_i_j; for i=1, 1 \le j \le nst
	
	JP_ij = J_i_i-1 Phi_i-1_j + J_i_i Phi_i_j, for 2 \le i \le Nt, 1 \le j \le nst
	
	Possible implementation: depending on the specific time integration scheme, (hence the sparsity of the huge space-time Jacobian), we can extract sequentially the specific spatial Jacobian at time steps n_t, n_t-1, etc. Note also that to compute spatial Jacobian at time step n_t, we'll need U(n_t) and/or U(n_t-1). 
	
	Data for this sub-step: U_st (uhat) to compute spatial jac, Phi_st, JP

	2.2.3 MultiSpaceTimeBasesMethods::MultiplySpaceTimeResidual
	Calculate RHS = JP^T * r_st (strictly for JP only)
	
	Size of JP is Ns x nst Nt, size of r_st is Ns x Nt, RHS is a small ROM vector with size nst x 1. This method creates the RHS vector on step 4 in section 2.3.1 page 3 in the report. 
	
	Data for this sub-step: U_st (uhat) to compute spatial residual at all time steps, JP	
	
	2.2.4 MultiSpaceTimeBasesMethods::MultiplyItself
	Calculate LHS = JP^T * JP (strictly for JP only)
	
	Size of JP is Ns x nst Nt, resulting LHS matrix is a small ROM matrix with size nst x nst. This method creates the LHS matrix on step 3 in section 2.3.1 page 3 in the report. 
	
	Data for this sub-step: JP
	*/

  // read from file the jacobian of the decoder
  constexpr int romSize = 1;
  // store modes computed before from file
  decoder_jac_t phi =
    pressio::apps::test::epetra::readBasis("basis.txt", romSize, 
    	numCell, Comm, appObj.getDataMap());

  if( phi.globalNumVectors() != romSize ) return 0;

  // create decoder obj
  decoder_t decoderObj(phi);

  // my reference state = a matrix Ns x Nt with all column vectors are initial (spatial) condition
  auto yRef = appObj.getInitialSpaceTimeState();

  // define ROM state and set to zero
  lspg_state_t yROM(romSize);
  yROM.putScalar(0.0);


	//====== step 3: create the space-time problem ======//
	/*--- step 3.1: here is where we instantiate an object of SpaceTimeProblem class 
	Define yROM, initial condition/reference solution	
	*/

  // define LSPG type
  using lspg_problem_type = pressio::rom::DefaultLSPGSpaceTimeTypeGenerator<
    fom_t, decoder_t, lspg_state_t>;

  pressio::rom::LSPGSpaceTimeProblemGenerator<lspg_problem_type> lspgProblem(
      appObj, *yRef, decoderObj, yROM);

  using rom_system_t = typename lspg_problem_type::lspg_system_t;


	//====== step 4: create a solver object ======//
	/* 
	declare/define linear solver to solve the resulting linear system
	declare nonlinear solver (e.g., Gauss-Newton) to solve optimization problem
	*/

  // linear solver
  using eig_dyn_mat  = Eigen::Matrix<scalar_t, -1, -1>;
  using hessian_t  = pressio::containers::Matrix<eig_dyn_mat>;
  using solver_tag   = pressio::solvers::linear::iterative::LSCG;
  using linear_solver_t = pressio::solvers::iterative::EigenIterative<solver_tag, hessian_t>;
  linear_solver_t linSolverObj;

  // (nonlinear) Gauss-Newton solver
  // hessian comes up in GN solver, it is (J phi)^T (J phi)
  // rom is solved using eigen, hessian is wrapper of eigen matrix
  using eig_dyn_mat  = Eigen::Matrix<scalar_t, -1, -1>;
  using converged_when_t = pressio::solvers::iterative::default_convergence;
  using gnsolver_t   = pressio::solvers::iterative::GaussNewton<
    rom_system_t, converged_when_t, linear_solver_t>;
  gnsolver_t solver(lspgProblem.systemObj_, yROM, linSolverObj);
  solver.setTolerance(1e-14);
  solver.setMaxIterations(200);


	//====== step 5: solve the space-time problem ======//
  solver.solve(lspgProblem.systemObj_, yROM);


	//====== step 6: check the solution of space-time problem ======//		
	/*
	Compute relative error between space-time solution with respect to FOM solution
	*/
  // reconstruct the fom corresponding to our rom final state (should be a matrix Ns x Nt)
  auto yFomApprox = lspgProblem.yFomReconstructor_(yROM);
  appObj.printStateToFile("rom.txt", *yFomApprox.data());
  {
  	// TODO: extract the last column of matrix yFomApprox and compare it with the FOM solution obtained with implicit euler, same time step, for 10 steps?

  }

  MPI_Finalize();
  std::cout << checkStr <<  std::endl;
  return 0;  
}

