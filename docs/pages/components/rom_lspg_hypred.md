
# rom: LSPG: hyper-reduced problem


Defined in: `<pressio/rom_lspg.hpp>`

Public namespace: `pressio::rom::lspg`

## Overview


## 1. Creating a problem instance

```cpp
// overload for continuous-time systems
template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,															  (1)
  class FomReferenceStateType,
  class HypRedOperatorUpdaterType
  >
ReturnType create_hyperreduced_unsteady_problem(pressio::ode::StepScheme,
												const FomSystemType &,
												DecoderType & ,
												const RomStateType &,
												const FomReferenceStateType &,
												const HypRedOperatorUpdaterType &)

// overload for discrete-time systems
template<
  std::size_t num_states,
  class FomSystemType,
  class DecoderType,
  class RomStateType,															  (2)
  class FomReferenceStateType
  >
ReturnType create_hyperreduced_unsteady_problem(const FomSystemType & fomSysObj,
												DecoderType & decoder,
												const RomStateType & stateIn,
												const FomReferenceStateType & fomRef);
```

### Parameters and Requirements

- `FomSystemType`:
  - your adapter class type specifying the FOM problem. <br/>
  - for 1: must satisfy the [continuous-time API](./md_pages_components_rom_fom_apis.html)
  - for 2: must satisfy the [discrete-time API](./md_pages_components_rom_fom_apis.html)

- `DecoderType`:
  - decoder class type
  - must satify the requirements listed [here](md_pages_components_rom_decoder.html)

- `RomStateType`:
  - currently, it must be either an Eigen vector or a Kokkos 1D view

- `FomReferenceStateType`:
  - your FOM reference state that is used when reconstructing the FOM state
  - must be copy-constructible and the following must be true:<br/>
  ```cpp
  std::is_same<FomReferenceStateType, typename DecoderType::fom_state_type>::value == true
  ```

- `HypRedOperatorUpdaterType`:
  - an instance of class that knows how to update operands that live on the stencil and sample mesh.
  - Must meet the following interface:

  ```cpp
  struct HypRedUpdater
  {
    using operand_type1 = /*should be the type of your fom velocity*/;
    using operand_type2 = /*should be the type of your decoder's jacobian*/;

	template<class ScalarType>
	void updateSampleMeshOperandWithStencilMeshOne(operand_type1 & operand_on_sample,
												   ScalarType alpha,
												   const operand_type1 & operand_on_stencil,
												   ScalarType beta) const
	{
	  // appropriately compute:
	  // operand_on_sample = alpha*operand_on_sample + beta*operand_on_stencil
	}

	template<class ScalarType>
	void updateSampleMeshOperandWithStencilMeshOne(operand_type2 & operand_on_sample,
												   ScalarType alpha,
												   const operand_type2 & operand_on_stencil,
												   ScalarType beta) const
	{
	  // appropriately compute:
	  // operand_on_sample = alpha*operand_on_sample + beta*operand_on_stencil
	}
  };
  ```

- `num_states`:
  - *total* number of states you need to use (must be <= 3), if you need more open issue
  - only needed for the discrete-time case


### Problem class API

A problem meets the following interface:

```cpp
class UnsteadyLspgProblem
{
public:
  using traits = /* nested typedef with trait class */;

  // returns the underlying stepper to use to solve the problem
  auto & stepper();

  // const ref to the object knowing how to reconstruct a FOM state
  const auto & fomStateReconstructor() const;
};
```
