
#ifndef PRESSIO_ROMS_TESTS_PREC_LSPG_HELPERS_HPP_
#define PRESSIO_ROMS_TESTS_PREC_LSPG_HELPERS_HPP_

struct MyCustomDecoder
{
  using jacobian_type  = pressio::containers::MultiVector<Eigen::MatrixXd>;
  using fom_state_type = ::pressio::containers::Vector<Eigen::VectorXd>;

private:
  mutable jacobian_type jac_;

public:
  MyCustomDecoder() = delete;

  MyCustomDecoder(const int fomSize, const int romSize)
    : jac_(fomSize, romSize)
  {
    jac_.data()->setConstant(1.);
  }

  const jacobian_type & jacobianCRef() const{
    return jac_;
  }

  template <typename rom_state_type>
  void updateJacobian(const rom_state_type &) const{}

  template <typename rom_state_type>
  void applyMapping(const rom_state_type & romState,
		  fom_state_type & result) const
  {
    result.data()->setConstant(1.);
  }
};

class Preconditioner
{
  using scalar_type = double;
  using state_type  = Eigen::VectorXd;
  using velocity_type = state_type;
  using dense_matrix_type = Eigen::MatrixXd;

public:
  // for steady
  void applyPreconditioner(const state_type & yState,
			   velocity_type & operand) const
  {
    for (auto i=0; i<operand.size(); ++i)
      operand(i)+=1.;
  }

  void applyPreconditioner(const state_type & yState,
			   dense_matrix_type & operand) const
  {
    assert(operand.cols()==3);
    for (auto i=0; i<operand.rows(); ++i)
      for (auto j=0; j<operand.cols(); ++j)
	operand(i,j) += static_cast<scalar_type>(i);
  }

  // for unsteady
  void applyPreconditioner(const state_type & yState,
			   scalar_type time,
			   velocity_type & operand) const
  {
    for (auto i=0; i<operand.size(); ++i)
      operand(i)+=1.;
  }

  void applyPreconditioner(const state_type & yState,
			   scalar_type time,
			   dense_matrix_type & operand) const
  {
    assert(operand.cols()==3);
    for (auto i=0; i<operand.rows(); ++i)
      for (auto j=0; j<operand.cols(); ++j)
	operand(i,j) += static_cast<scalar_type>(i);
  }

};

#endif
