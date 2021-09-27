
# rom: LSPG: masked problem


@m_class{m-note m-default}

@parblock
Defined in: `<pressio/rom_lspg.hpp>`

Public namespace: `pressio::rom::lspg`
@endparblock


## Overview


## 1. Creating a problem instance

```cpp
// overload for continuous-time systems
template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,															  (1)
  class FomReferenceStateType,
  class MaskerType
  >
ReturnType create_masked_unsteady_problem(pressio::ode::StepScheme,
										  const FomSystemType &,
										  DecoderType &,
										  const RomStateType &,
										  const FomReferenceStateType &,
										  const MaskerType &);

// overload for discrete-time systems
template<
  std::size_t num_states,
  class FomSystemType,
  class DecoderType,
  class RomStateType,															  (2)
  class FomReferenceStateType,
  class MaskerType
  >
ReturnType create_masked_unsteady_problem(const FomSystemType &,
										  DecoderType &,
										  const RomStateType &,
										  const FomReferenceStateType &,
										  const MaskerType &);
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

- `RomStateType`:
  - currently, it must be either an Eigen vector or a Kokkos 1D view

- `num_states`:
  - *total* number of states you need to use (must be <= 3), if you need more open issue
  - only needed for the discrete-time case

- `MaskedType`:
  - an functor responsible of "masking" the FOM operators
  - must be a functor with a specific API:

```cpp
struct ValidMasker
{
  using operand_type1 = /*same type as your fom velocity*/;
  using operand_type2 = /*same type as your decoder's jacobian*/;

  operand_type1 createApplyMaskResult(const operand_type1 & unmasked_object);
  operand_type2 createApplyMaskResult(const operand_type2 & unmasked_object);

  template<class TimeType>
  void operator()(const operand_type1 & unmasked_object,
                  const TimeType time,
				  operand_type1 & result);

  template<class TimeType>
  void operator()(const operand_type2 & unmasked_object,
                  const TimeType time,
				  operand_type2 & result);
};
```



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
