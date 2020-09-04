
#ifndef PRESSIO_TESTS_BURG1D_CUSTOM_MAPPING_HPP_
#define PRESSIO_TESTS_BURG1D_CUSTOM_MAPPING_HPP_

#include "pressio_rom.hpp"
#include "utils_eigen.hpp"

struct MyCustomDecoder
{
  // this is mandatory because pressio detects it
  using jacobian_type  = pressio::containers::MultiVector<Eigen::MatrixXd>;

private:
  const int romSize_ = {};
  mutable jacobian_type jac_;
  mutable int stepsCount1_ = 0;
  mutable int stepsCount2_ = 0;

public:
  MyCustomDecoder() = delete;

  MyCustomDecoder(const int fomSize, const int romSize)
    : romSize_{romSize}, jac_(fomSize, romSize)
  {
    //initialize the jacobian to be [1,2,3]
    for (int i=0; i<fomSize; ++i){
      jac_(i,0) = 1; jac_(i,1) = 2; jac_(i,2) = 3;
      jac_(i,3) = 4; jac_(i,4) = 5;
    }
  }

  template <typename rom_state_type>
  void updateJacobian(const rom_state_type &) const
  {
    ++stepsCount2_;
    // here, we should have the counts to be equal if things are right
    if (stepsCount1_ != stepsCount2_){
      throw std::runtime_error
	("invalid logic for Galerkin calls tocustomDecoder, something is off");
    }

    // at step 2, change the jacobian
    if(stepsCount2_ == 2){
      for (int i=0; i<jac_.extent(0); ++i){
    	jac_(i,0) = 2; jac_(i,1) = 3; jac_(i,2) = 4;
    	jac_(i,3) = 5; jac_(i,4) = 6;
      }
    }
    // at step 3, change the jacobian again
    if(stepsCount2_ == 3){
      for (int i=0; i<jac_.extent(0); ++i){
    	jac_(i,0) = 3; jac_(i,1) = 4; jac_(i,2) = 5;
    	jac_(i,3) = 6; jac_(i,4) = 7;
      }
    }
  }

  template <typename rom_state_type>
  void applyMapping(const rom_state_type & romState,
		    ::pressio::containers::Vector<Eigen::VectorXd> & result) const
  {
    ++stepsCount1_;
    // applyMapping should be called BEFORE the updateJacobian
    // because pressio first reconstrct FOM, then call FOM velocity
    // and then project the velocity. // Here we can use counters to verify this
    if (stepsCount1_ != stepsCount2_+1){
      throw std::runtime_error
	("invalid logic for Galerkin calls tocustomDecoder, something is off");
    }

    const auto & jacNativeObj = *jac_.data();
    const auto & romStateNativeObj = *romState.data();
    auto & resultNativeObj = *result.data();
    resultNativeObj = jacNativeObj * romStateNativeObj;
  }

  const jacobian_type & getReferenceToJacobian() const{
    return jac_;
  }


};//end

#endif
