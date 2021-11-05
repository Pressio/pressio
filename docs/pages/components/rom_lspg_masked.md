
# rom: LSPG: masked problem


@m_class{m-note m-default}

@parblock
Defined in: `<pressio/rom_lspg.hpp>`

Public namespace: `pressio::rom::lspg`
@endparblock

<br/>

@m_class{m-block m-warning}

@par Prerequisite reading:
Before you read this page, make sure you
read the [overall design idea of the unsteady LSPG](md_pages_components_rom_lspg_unsteady.html).
@endparblock

## API

```cpp
// overload for continuous-time systems
template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,															  (1)
  class FomReferenceStateType,
  class MaskerType
  >
ReturnType create_masked_unsteady_problem(pressio::ode::StepScheme scheme,
										  const FomSystemType & fomSystem,
										  DecoderType & decoder,
										  const RomStateType & romState,
										  const FomReferenceStateType & fomRefState,
										  const MaskerType & masker);

// overload for discrete-time systems
template<
  std::size_t num_states,
  class FomSystemType,
  class DecoderType,
  class RomStateType,															  (2)
  class FomReferenceStateType,
  class MaskerType
  >
ReturnType create_masked_unsteady_problem(const FomSystemType & fomSystem,
										  DecoderType & decoder,
										  const RomStateType & romState,
										  const FomReferenceStateType & fomRefState,
										  const MaskerType & masker);
```

### Parameters and Requirements

- `fomSystem`:
  - instance of your adapter class type specifying the FOM problem
  - for 1: must satisfy the [continuous-time API](./md_pages_components_rom_fom_apis.html)
  - for 2: must satisfy the [discrete-time API](./md_pages_components_rom_fom_apis.html)

- `decoder`:
  - decoder object
  - must satify the requirements listed [here](md_pages_components_rom_decoder.html)

- `romState`:
  - currently, it must be either an Eigen vector or a Kokkos 1D view

- `fomRefState`:
  - your FOM reference state that is used when reconstructing the FOM state
  - must be copy-constructible and the following must be true:<br/>
  ```cpp
  std::is_same<FomReferenceStateType, typename DecoderType::fom_state_type>::value == true
  ```

- `romState`:
  - currently, it must be either an Eigen vector or a Kokkos 1D view

- `num_states`:
  - *total* number of states you need to use (must be <= 3), if you need more open issue
  - only needed for the discrete-time case

- `masker`:
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
