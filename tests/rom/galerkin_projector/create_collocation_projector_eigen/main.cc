
#include "pressio_rom_galerkin.hpp"

using eigen_v_t = Eigen::VectorXd;
using eigen_m_t = Eigen::MatrixXd;
using fom_state_t	 = pressio::containers::Vector<eigen_v_t>;
using decoder_jacobian_t = pressio::containers::MultiVector<eigen_m_t>;
using decoder_t		 = pressio::rom::LinearDecoder<decoder_jacobian_t, fom_state_t>;

struct MyCollocator
{
  std::vector<int> rows_;

  MyCollocator(std::initializer_list<int> l) : rows_(l){}

  eigen_m_t sampleRows(const eigen_m_t & operand)
  {
    eigen_m_t result(rows_.size(), operand.cols());
    for (auto i=0; i<(int)rows_.size(); ++i){
      for (auto j=0; j<result.cols(); ++j){
	result(i,j) = operand(rows_[i], j);
      }
    }
    return result;
  }
};

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
  decoder_jacobian_t phi(10,3);
  for (auto i=0; i<phi.extent(0); ++i)
    for (auto j=0; j<phi.extent(1); ++j)
      phi(i,j) = (double) i;

  // create decoder
  decoder_t decoder(phi);

  // only pick rows 0,3,4,5,8
  MyCollocator mapper({0,3,4,5,8});
  /*
    projector represents this operator
    matrix = [0 0 0;
	      3 3 3;
              4 4 4;
	      5 5 5;
	      8 8 8]
  */
  using pressio::rom::galerkin::createCollocationProjector;
  auto galProjector = createCollocationProjector(decoder, std::move(mapper));

  // create a vector to test project
  // vector must have 5 entries, i.e. same as collocation points
  fom_state_t operand(5);
  pressio::ops::fill(operand, 2.);

  using galerkin_rhs_t = pressio::containers::Vector<eigen_v_t>;
  galerkin_rhs_t r(3);

  /* compute matrix^T operand*/
  galProjector.apply(operand, r);

  /* result should be [40 40 40] */
  const std::vector<double> gold = {40., 40., 40.};

  // check the reconstructed fom state
  for (auto i=0; i<r.extent(0); i++){
    if (std::abs(gold[i] - r(i)) > 1e-13){
      checkStr = "FAILED";
      break;
    }
  }

  std::cout << checkStr <<  std::endl;
  pressio::log::finalize();
  return 0;
}
