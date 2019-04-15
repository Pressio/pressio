#include "CORE_ALL"
#include "APPS_STEADYLINADVDIFF2D"
#include "../fom_gold_states.hpp"

constexpr double eps = 1e-7;
std::string checkStr {"PASSED"};

template <typename T>
void checkSol(const T & y,
	      const std::map<int, double> & trueS){
  const auto map = y.Map();
  auto * globID = map.MyGlobalElements();
  const auto myN = map.NumMyElements();
  for(auto i=0; i<myN; i++){
    if (std::abs(y[i] - trueS.at(globID[i])) > eps) checkStr = "FAILED";
  }
}

int main(int argc, char *argv[]){
  using fom_t	 = ::rompp::apps::SteadyLinAdvDiff2dEpetra;

  int rank;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);

  const int Nx = 11, Ny = Nx*2-1;
  fom_t  appObj(Comm, Nx, Ny);
  appObj.setup();
  appObj.assembleMatrix();
  appObj.fillRhs();
  appObj.solve();
  //appObj.printStateToFile("out.txt");
  auto T = appObj.getState();
  T->Print(std::cout << std::setprecision(15));
  checkSol(*T, rompp::apps::test::steadyAdvDiff2d_nx11ny21);

  MPI_Finalize();
  std::cout << checkStr <<  std::endl;
  return 0;
}
