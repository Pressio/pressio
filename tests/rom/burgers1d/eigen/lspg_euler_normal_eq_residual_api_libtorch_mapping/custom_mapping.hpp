
#ifndef PRESSIO_TESTS_BURG1D_CUSTOM_MAPPING_HPP_
#define PRESSIO_TESTS_BURG1D_CUSTOM_MAPPING_HPP_

#include "pressio_rom.hpp"
#include "utils_eigen.hpp"
#include <torch/script.h>

template <
  typename native_dense_mat_type, // this is type of native dense matrix
  typename fom_state_type  // this is a pressio::containers::Vector<>
  >
struct MyCustomDecoder
{
  // this is mandatory because pressio detects it
  using jacobian_type  = pressio::containers::MultiVector<native_dense_mat_type>;

private:
  const int romSize_ = {};
  mutable jacobian_type jac_;

public:
  MyCustomDecoder() = delete;

  MyCustomDecoder(const int romSize, const int numCell)
    : romSize_{romSize},
      jac_{numCell, romSize}
      //jac_{pressio::rom::test::eigen::readBasis("basis.txt", romSize, numCell)}
  {}

  template <typename rom_state_type>
  void applyMapping(const rom_state_type & romState, fom_state_type & result) const
  {
    // here romState has same type as the one you used for lspg_state_t in main
    // result is a pressio::containers::Vector< native_state_type>
    // to get a reference to native data, use data() as follows:

    // Create module here because applyMapping() is required to be const and module.forward is nonstatic 
    torch::jit::script::Module module;
    try {
        // Deserialize the ScriptModule from a file using torch::jit::load().
        // TODO: get rid of hardcoded path
        module = torch::jit::load("/scratch/cpp_autoencoder/burgers1d/traced_model.pt");
        }
    catch (const c10::Error& e) {
        std::cerr << "error loading the model\n";
        }
    
    // Create a vector of inputs.
    std::vector<torch::jit::IValue> inputs;
    torch::Tensor f=torch::empty({1,romSize_});
    for(int i=0; i<romSize_; i++)
      f[0][i] = float(romState[i]);
    inputs.push_back(f);

    // Execute model
    at::Tensor model_output = module.forward(inputs).toTensor();
    int numCell = model_output.sizes()[2];

    // Save output
    for(int i=0; i<numCell; i++)
      result[i] = model_output[0][0][i].item<double>();

    // Update Jacobian
    updateJacobian(romState);
  }

  template <typename rom_state_type>
  void updateJacobian(const rom_state_type & romState) const
  {
    torch::jit::script::Module module;
    try {
        // Deserialize the ScriptModule from a file using torch::jit::load().
        // TODO: get rid of hardcoded path
        module = torch::jit::load("/scratch/cpp_autoencoder/burgers1d/traced_model.pt");
        }
    catch (const c10::Error& e) {
        std::cerr << "error loading the model\n";
        }

    // Create a vector of inputs.
    std::vector<torch::jit::IValue> inputs;
    torch::Tensor f=torch::empty({1,romSize_});
    for(int i=0; i<romSize_; i++)
      f[0][i] = float(romState[i]);
    inputs.push_back(f);

    // Execute model for romState
    at::Tensor model_output = module.forward(inputs).toTensor();
    int numCell = model_output.sizes()[2];

    // Execute model for perturbed romStates and update Jacobian
    const float eps = 0.01;
    for(int param=0; param<romSize_; param++) {
      f[0][param] += eps;
      at::Tensor shifted_model_output = module.forward(inputs).toTensor();
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
