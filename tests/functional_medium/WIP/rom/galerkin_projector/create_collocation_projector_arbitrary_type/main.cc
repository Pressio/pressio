
#include "pressio_rom_galerkin.hpp"

using scalar_t = double;
using v_t = std::vector<scalar_t>;
using m_t = std::vector<v_t>;
using fom_state_t	 = pressio::containers::Vector<v_t>;
using decoder_jacobian_t = pressio::containers::DenseMatrix<m_t>;
using decoder_t		 = pressio::rom::LinearDecoder<decoder_jacobian_t, fom_state_t>;

struct myOps
{
  //  y = beta * y + alpha*A^T*x
  template <typename y_t>
  void product(::pressio::transpose mode,
	       scalar_t alpha,
	       const m_t & A,
	       const v_t & x,
	       scalar_t beta,
	       y_t & y) const
  {
    // y is sharedmem 1d array and subscriptable with (i)
    const auto nArows = A.size();
    const auto nAcols = A.cbegin()->size();
    for (std::size_t j=0; j<nAcols; ++j){
      y(j) = beta*y(j);
      for (std::size_t i=0; i<nArows; ++i)
	y(j) += alpha * A[i][j] * x[i];
    }
  }

  //  C = beta * C + alpha*A^T*B
  template <typename C_t>
  void product(::pressio::transpose,
	       ::pressio::nontranspose,
	       const scalar_t alpha,
	       const m_t & A,
	       const m_t & B,
	       const scalar_t beta,
	       C_t & C) const
  {
    // C is sharedmem matrix and subscriptable with (i,j)
    const auto nArows = A.size();
    const auto nAcols = A.cbegin()->size();
    const auto nCrows = C.extent(0);
    const auto nCcols = C.extent(1);
    for (std::size_t i=0; i<nArows; ++i){
      for (std::size_t j=0; j<nAcols; ++j)
      {
	C(i,j) = beta*C(i,j);
	scalar_t sum = {};
	for (std::size_t k=0; k<nArows; ++k){
	  sum += A[k][j] * B[k][j];
	}
	C(i,j) += alpha*sum;
      }
    }
  }
};

struct MyCollocator
{
  std::vector<int> rows_;

  MyCollocator(std::initializer_list<int> l) : rows_(l){}

  m_t sampleRows(const m_t & operand)
  {
    const auto nCols = operand.cbegin()->size();

    m_t result(rows_.size());
    for (auto & it : result) it.resize(nCols);

    const auto numCols = operand.cbegin()->size();
    for (auto i=0; i<(int)rows_.size(); ++i)
      for (auto j=0; j<(int)numCols; ++j)
	result[i][j] = operand[rows_[i]][j];
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
  m_t phi0(10); for (auto & it : phi0) it.resize(3);
  for (auto i=0; i<(int)phi0.size(); ++i)
    for (auto j=0; j<(int)phi0[i].size(); ++j)
      phi0[i][j] = (scalar_t) i;

  decoder_jacobian_t phi(phi0);
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
  myOps ops;
  using pressio::rom::galerkin::createCollocationProjector;
  auto galProjector = createCollocationProjector(decoder, std::move(mapper), ops);

  // create a vector to test project
  // vector must have 5 entries, i.e. same as collocation points
  fom_state_t operand(v_t(5, 2.));

  using galerkin_rhs_t = pressio::containers::Vector<Eigen::VectorXd>;
  galerkin_rhs_t r(3);

  /* compute matrix^T operand*/
  galProjector.apply(operand, r);

  /* result should be [40 40 40] */
  const std::vector<scalar_t> gold = {40., 40., 40.};

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
