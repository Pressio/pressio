#include "CONTAINERS_ALL"
#include "APPS_STEADYLINADVDIFF2D"
#include <random>
#include <Eigen/SVD>

using gen_t	   = std::mt19937;
using rand_distr_t = std::uniform_real_distribution<double>;

constexpr double eps = 1e-7;
std::string checkStr {"PASSED"};

// range of Prandtl number
constexpr std::array<double,2> Pr_range{{1.0, 5.0}};

// range of Reynolds number
constexpr std::array<double,2> Re_range{{10., 100.0}};

void readMatrixFromFile(std::string filename,
			std::vector<std::vector<double>> & A0,
			int ncols){
  assert( A0.empty() );
  std::ifstream source;
  source.open( filename, std::ios_base::in);
  std::string line, colv;
  std::vector<double> tmpv(ncols);
  while (std::getline(source, line) ){
    std::istringstream in(line);
    for (int i=0; i<ncols; i++){
      in >> colv;
      tmpv[i] = atof(colv.c_str());
    }
    A0.emplace_back(tmpv);
  }
  source.close();
}

int main(int argc, char *argv[]){
  using fom_t	 = ::pressio::apps::SteadyLinAdvDiff2dEpetra;

  int rank;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  assert(Comm.NumProc() == 1);

  //---------------------------------------
  // generate random samples of Parameters

  // random number generator (seeded)
  unsigned int seed = 1343234343;
  std::mt19937 engine(seed);
  rand_distr_t distr(0., 1.0);
  auto genPr = [&distr, &engine](){
		 auto c1 = Pr_range[1]-Pr_range[0];
		 auto c2 = Pr_range[0];
		 return c1 * distr(engine) + c2;
	     };

  auto genRe = [&distr, &engine](){
		 auto c1 = Re_range[1]-Re_range[0];
		 auto c2 = Re_range[0];
		 return c1 * distr(engine) + c2;
	     };

  // number of sample to take
  constexpr int nSamples = 5;

  // fill for Prandtl
  std::vector<double> PrS(nSamples);
  std::generate(PrS.begin(), PrS.end(), genPr);

  // fill for Reynolds
  std::vector<double> ReS(nSamples);
  std::generate(ReS.begin(), ReS.end(), genRe);

  if(rank==0){
    auto it1 = PrS.begin();
    auto it2 = ReS.begin();
    for( ;it2<ReS.end(); it1++, it2++)
      std::cout << std::setprecision(15)
		<< *it1 << " " << *it2
		<< "\n";
  }

  //---------------------------------------
  // fix discretization for all samples
  const int Nx = 11, Ny = Nx*2-1;

  /* # of dofs is != Nx*Ny because of how we solve pdd */
  const int numDof = (Nx-2)*Ny;

  // create as many app objects as samples
  std::vector<fom_t> vecObjs;
  for (auto i=0; i<nSamples; i++){
    vecObjs.emplace_back(Comm,Nx, Ny, PrS[i], ReS[i]);
  }

  // solve all problems
  for (auto & it : vecObjs){
    it.setup();
    it.assembleMatrix();
    it.fillRhs();
    it.solve();
  }

  // collect all solutions into matrix
  // I can do this this easily because we know # ranks = 1
  using eig_mat = Eigen::MatrixXd;
  eig_mat A(numDof, nSamples);
  int j=0;
  for (const auto & it : vecObjs){
    auto T = it.getState();
    for (auto i=0; i<numDof; i++)
      A(i,j) = (*T)[i];
    j++;
  }

  // do SVD
  Eigen::JacobiSVD<eig_mat> svd(A, Eigen::ComputeThinU);
  auto U = svd.matrixU();
  std::cout << std::setprecision(15) << U << std::endl;

  // read gold basis from file
  std::vector<std::vector<double>> goldU;
  readMatrixFromFile("gold_basis.txt", goldU, nSamples);

  // check that computed matches gold
  assert( (size_t) goldU.size() == (size_t) U.rows() );
  for (auto i=0; i<U.rows(); i++)
    for (j=0; j<nSamples; j++)
      if ( std::abs(goldU[i][j] - U(i,j)) > eps ) checkStr = "FAILED";

  MPI_Finalize();
  std::cout << checkStr <<  std::endl;
  
  return 0;
}
