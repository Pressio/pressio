
# rom: LSPG: steady hyper-reduced problem


@m_class{m-note m-default}

@parblock
Defined in: `<pressio/rom_lspg.hpp>`

Public namespace: `pressio::rom::lspg`
@endparblock


<br/>

@m_class{m-block m-warning}

@par Prerequisite reading:
Before you read this page, make sure you
read the [overall design idea behind steady LSPG](md_pages_components_rom_lspg_steady.html).
@endparblock


## Creating a problem instance


```cpp
template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType
  >
auto create_hyperreduced_steady_problem(const FomSystemType & fomSystem,
										DecoderType & decoder,
										const RomStateType & romState,
										const FomReferenceStateType & fomRefState);

template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class PreconditionerType
  >
auto create_hyperreduced_steady_problem(const FomSystemType & fomSystem,
										DecoderType & decoder,
										const RomStateType & romState,
										const FomReferenceStateType & fomRefState,
										const PreconditionerType & preconditioner);
```

### Parameters and Requirements

- `fomSystem`:
  - instance of your FOM adapter specifying the FOM problem
  - must satisfy the steady API, see [here](./md_pages_components_rom_fom_apis.html)

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
