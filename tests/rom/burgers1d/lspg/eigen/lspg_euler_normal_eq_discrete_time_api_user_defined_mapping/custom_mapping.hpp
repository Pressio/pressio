
#ifndef PRESSIO_TESTS_BURG1D_CUSTOM_MAPPING_HPP_
#define PRESSIO_TESTS_BURG1D_CUSTOM_MAPPING_HPP_

#include "utils_eigen.hpp"

template <
  typename native_dense_mat_type,// this is type of native dense matrix
  typename fom_state_t// this is a pressio::containers::Vector<>
  >
struct MyCustomDecoder
{
  // this is mandatory because pressio detects it
  using jacobian_type = pressio::containers::MultiVector<native_dense_mat_type>;
  using fom_state_type = fom_state_t;

private:
  const int romSize_ = {};
  jacobian_type jac_;

public:
  MyCustomDecoder() = delete;

  MyCustomDecoder(const int romSize, const int numCell)
    : romSize_{romSize},
      jac_{pressio::rom::test::eigen::readBasis("basis.txt", romSize, numCell)}
  {
  }

  template <typename rom_state_type>
  void applyMapping(const rom_state_type & romState, fom_state_type & result) const
  {
    // here romState has same type as the one you used for lspg_state_t in main
    // result is a pressio::containers::Vector< native_state_type>
    // to get a reference to native data, use data() as follows:
    const auto & jacNativeObj = *jac_.data();
    const auto & romStateNativeObj = *romState.data();
    auto & resultNativeObj = *result.data();

    // having the native objects, you can do operation using native functionalities
    resultNativeObj = jacNativeObj * romStateNativeObj;
  }

  const jacobian_type & jacobianCRef() const
  {
    return jac_;
  }

  template <typename rom_state_type>
  void updateJacobian(const rom_state_type &)
  {
    // no op, the Jacobian is fixed
  }
};//end

#endif
