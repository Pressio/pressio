
#include "pressio_rom_galerkin.hpp"

using eigen_v_t = Eigen::VectorXd;
using eigen_m_t = Eigen::MatrixXd;
using fom_state_t   = pressio::containers::Vector<eigen_v_t>;
using proj_matrix_t = pressio::containers::MultiVector<eigen_m_t>;

int main(int argc, char *argv[])
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});
  std::string checkStr {"PASSED"};

  /*
     phi = [0 0 0;
	    1 1 1;
	    2 2 2;
	    ...
	    9 9 9]
  */
  proj_matrix_t phi(10,3);
  for (auto i=0; i<phi.extent(0); ++i)
    for (auto j=0; j<phi.extent(1); ++j)
      phi(i,j) = (double) i;

  using pressio::rom::galerkin::createArbitraryProjector;
  auto galProjector = createArbitraryProjector(std::move(phi));

  // create a vector to test project
  // vector must have 5 entries, i.e. same as collocation points
  fom_state_t operand(10);
  pressio::ops::fill(operand, 1.);

  using galerkin_rhs_t = pressio::containers::Vector<eigen_v_t>;
  galerkin_rhs_t r(3);

  /* compute matrix^T operand*/
  galProjector.apply(operand, r);

  /* result should be [45 45 45] */
  std::vector<double> gold(10, 45.);

  // check the reconstructed fom state
  for (auto i=0; i<r.extent(0); i++){
    if (std::abs(gold[i] - r(i)) > 1e-13){
      checkStr = "FAILED";
      break;
    }
  }

  std::cout << checkStr <<  std::endl;
  return 0;
}
