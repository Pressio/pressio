
# rom: LSPG: steady masked problem


@m_class{m-note m-default}

@parblock
Defined in: `<pressio/rom_lspg.hpp>`

Public namespace: `pressio::rom::lspg`
@endparblock



## Overview
\todo


## 1. Creating a problem instance

```cpp
template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class MaskerType
  >
auto create_masked_steady_problem(const FomSystemType &,
								  DecoderType &,
								  const RomStateType &,
								  const FomReferenceStateType &,
								  const MaskerType &);

template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class MaskerType,
  class PreconditionerType
  >
auto create_masked_steady_problem(const FomSystemType &,
								  DecoderType &,
								  const RomStateType &,
								  const FomReferenceStateType &,
								  const PreconditionerType &,
								  const MaskerType &);
```

### Parameters and Requirements

- `FomSystemType`:
  - your adapter class type specifying the FOM problem. <br/>
  - Must satisfy the steady API, see [here](./md_pages_components_rom_fom_apis.html)

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


### Problem class API

A problem meets the following interface:

```cpp
class SteadyLspgProblem
{
public:
  using traits = /* nested typedef with trait class */;

  // returns the underlying system to use to solve the problem
  auto & system();

  // const ref to the object knowing how to reconstruct a FOM state
  const auto & fomStateReconstructor() const;
};
```
