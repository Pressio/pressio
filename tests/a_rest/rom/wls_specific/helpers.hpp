
#ifndef PRESSIO_ROMS_TESTS_WLS_CHECK_NORMALEQNS_HELPERS_HPP_
#define PRESSIO_ROMS_TESTS_WLS_CHECK_NORMALEQNS_HELPERS_HPP_

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
  void updateJacobian(const rom_state_type &){}

  template <typename rom_state_type>
  void applyMapping(const rom_state_type & romState,
		  fom_state_type & result) const
  {
    result.data()->setConstant(1.);
  }
};

#endif
