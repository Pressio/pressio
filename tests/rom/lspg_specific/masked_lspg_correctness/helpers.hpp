
#ifndef PRESSIO_ROMS_TESTS_PREC_LSPG_HELPERS_HPP_
#define PRESSIO_ROMS_TESTS_PREC_LSPG_HELPERS_HPP_

struct MyCustomDecoder
{
  using jacobian_type = pressio::containers::MultiVector<Eigen::MatrixXd>;
  using fom_state_type = ::pressio::containers::Vector<Eigen::VectorXd>;

private:
  jacobian_type jac_;

public:
  MyCustomDecoder() = delete;

  MyCustomDecoder(const int fomSize, const int romSize)
    : jac_(fomSize, romSize)
  {
    jac_.data()->setConstant(1.);
  }

  const jacobian_type & jacobianCRef() const
  {
    return jac_;
  }

  template <typename rom_state_type>
  void updateJacobian(const rom_state_type &) {}

  template <typename rom_state_type>
  void applyMapping(const rom_state_type & romState,
		    fom_state_type & result) const
  {
    *result.data() = (*jac_.data()) * (*romState.data());

    std::cout << "applyMap result " << std::endl;
    std::cout << *result.data() << std::endl;
  }
};

class Masker
{
  int maskSize_ = {};
  std::vector<int> indices_ = {};

  using scalar_type = double;
  using state_type = Eigen::VectorXd;
  using residual_type = state_type;
  using velocity_type = state_type;
  using dense_matrix_type = Eigen::MatrixXd;

public:
  Masker(std::vector<int> indices)
    : maskSize_(indices.size()),
      indices_(indices)
  {
    assert(maskSize_ == 5);
  }

  residual_type createApplyMaskResult(const residual_type & operand) const
  {
    return residual_type(maskSize_);
  }

  dense_matrix_type createApplyMaskResult(const dense_matrix_type & operand) const
  {
    return dense_matrix_type(maskSize_, operand.cols());
  }

  // for steady
  void applyMask(const residual_type & operand,
		 residual_type & result) const
  {
    for(auto i = 0; i < result.size(); ++i)
      result(i) = operand(indices_[i]);
  }

  void applyMask(const dense_matrix_type & operand,
		 dense_matrix_type & result) const
  {
    for(auto i = 0; i < result.rows(); ++i)
      for(auto j = 0; j < result.cols(); ++j)
	result(i, j) = operand(indices_[i], j);
  }

  // for unsteady
  void applyMask(const residual_type & operand,
		 scalar_type time,
		 residual_type & result) const
  {
    for(auto i = 0; i < result.size(); ++i)
      result(i) = operand(indices_[i]);
  }

  void applyMask(const dense_matrix_type & operand,
		 scalar_type time,
		 dense_matrix_type & result) const
  {
    for(auto i = 0; i < result.rows(); ++i)
      for(auto j = 0; j < result.cols(); ++j)
	result(i, j) = operand(indices_[i], j);
  }
};

#endif
