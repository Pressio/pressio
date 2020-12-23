
#ifndef PRESSIO_TESTS_BURG1D_CUSTOM_MAPPING_HPP_
#define PRESSIO_TESTS_BURG1D_CUSTOM_MAPPING_HPP_

// #include "pressio_rom.hpp"
#include "utils_eigen.hpp"
#include <torch/script.h>

template <
  typename native_dense_mat_type, // this is type of native dense matrix
  typename fom_state_t  // this is a pressio::containers::Vector<>
  >
struct MyCustomDecoder
{
  // this is mandatory because pressio detects it
  using jacobian_type  = pressio::containers::MultiVector<native_dense_mat_type>;
  using fom_state_type = fom_state_t;

private:
  const int romSize_ = {};
  mutable jacobian_type jac_;
  const char* filename_;
  mutable torch::jit::script::Module decoder_;
public:
  MyCustomDecoder() = delete;

  MyCustomDecoder(const int romSize, const int numCell, const char* filename)
    : romSize_{romSize},
      jac_{numCell, romSize},
      filename_{filename}
  {
    try {
        decoder_ = torch::jit::load(filename);
        }
    catch (const c10::Error& e) {
        std::cerr << "error loading the model\n";
        }

  }

  template <typename rom_state_type>
  void applyMapping(const rom_state_type & romState, fom_state_type & result) const
  {
    // here romState has same type as the one you used for lspg_state_t in main
    // result is a pressio::containers::Vector< native_state_type>

    // Create a vector of inputs.
    std::vector<torch::jit::IValue> inputs;
    torch::Tensor f=torch::empty({1,romSize_});
    for(int i=0; i<romSize_; i++)
      f[0][i] = float(romState[i]);
    inputs.push_back(f);

    // Execute model
    at::Tensor model_output = decoder_.forward(inputs).toTensor();
    int numCell = model_output.sizes()[2];

    // Save output
    for(int i=0; i<numCell; i++)
      result[i] = model_output[0][0][i].item<double>();

    // Update Jacobian
    // TODO: move to getReferenceToJacobian() once romState is accessible there
    updateJacobian(romState);
  }

  template <typename rom_state_type>
  void updateJacobian(const rom_state_type & romState) const
  {
    // Create a vector of inputs.
    std::vector<torch::jit::IValue> inputs;
    torch::Tensor f=torch::empty({1,romSize_});
    for(int i=0; i<romSize_; i++)
      f[0][i] = float(romState[i]);
    inputs.push_back(f);

    // Execute model for romState
    at::Tensor model_output = decoder_.forward(inputs).toTensor();
    int numCell = model_output.sizes()[2];

    // Execute model for perturbed romStates and update Jacobian
    const float eps = 0.01;
    for(int param=0; param<romSize_; param++) {
      f[0][param] += eps;
      at::Tensor shifted_model_output = decoder_.forward(inputs).toTensor();
      f[0][param] -= eps;

      // Update Jacobian column
      for(int i=0; i<numCell; i++)
        jac_(i,param) = ( shifted_model_output[0][0][i].item<double>() - model_output[0][0][i].item<double>() )/double(eps);
    }
  }

  const jacobian_type & getReferenceToJacobian() const{
    return jac_;
  }
};//end

#endif
