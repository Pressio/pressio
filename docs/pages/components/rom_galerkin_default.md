
# rom: Galerkin: default problem


@m_class{m-note m-default}

@parblock
Defined in: `<pressio/rom_galerkin.hpp>`

Public namespace: `pressio::rom::galerkin`
@endparblock

<br/>

@m_class{m-block m-warning}

@par Prerequisite reading:
Before you read this page, make sure you
read the [overall design idea of Galerkin](md_pages_components_rom_galerkin.html).
@endparblock

## API

```cpp
template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType
  >																				 (1)
auto create_default_explicit_problem(pressio::ode::StepScheme scheme,
								     const FomSystemType & fomSystem,
									 DecoderType & decoder,
									 const RomStateType & romState,
									 const FomReferenceStateType & fomRefState)

template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType
  >																				 (2)
auto create_default_implicit_problem(pressio::ode::StepScheme scheme,
								     const FomSystemType & fomSystem,
									 DecoderType & decoder,
									 const RomStateType & romState,
									 const FomReferenceStateType & fomRefState)
```

This function returns an instance of the desired Galerkin problem.

### Parameters and Requirements

- `scheme`:
  - enum value to specify the stepper scheme
  - (1) explicit Galerkin, see [valid enum scheme choices](md_pages_components_ode_steppers_explicit.html)
  - (2) implicit Galerkin, seee [valid enum scheme choices](md_pages_components_ode_steppers_implicit.html)

- `fomSystem`:
  - instance of your FOM adapter specifying the FOM problem
  - Must satisfy one of the APIs suitable for Galerkin, see [API list](./md_pages_components_rom_fom_apis.html)

- `decoder`:
  - your decoder object
  - must satify the requirements listed [here](md_pages_components_rom_decoder.html)

- `romState`:
  - currently, it must be either an Eigen vector or a Kokkos 1D view

- `fomRefState`:
  - your FOM reference state that is used when reconstructing the FOM state
  - must be copy-constructible and the following must be true:<br/>
  ```cpp
  std::is_same<FomReferenceStateType, typename DecoderType::fom_state_type>::value == true
  ```

<br/>

@m_div{m-button m-primary}
<a href="https://pressio.github.io/pressio-tutorials/known_types/proms_galerkin_default_explicit.html">
@m_div{m-big}Explicit Galerkin sample code 1: @m_enddiv
@m_div{m-small} (for data types natively supported) @m_enddiv
</a> @m_enddiv

<!-- @m_div{m-button m-primary}
<a href="https://pressio.github.io/pressio-tutorials/html/md_pages_rom_galerkin_default_explicit_custom_types.html">
@m_div{m-big} Explicit Galerkin sample code 2: @m_enddiv
@m_div{m-small} (for arbitrary data types) @m_enddiv
</a> @m_enddiv
 -->
@m_div{m-button m-primary}
<a href="https://pressio.github.io/pressio-tutorials/known_types/proms_galerkin_default_implicit.html">
@m_div{m-big}Implicit Galerkin sample code 1: @m_enddiv
@m_div{m-small} (for data types natively supported) @m_enddiv
</a> @m_enddiv

<!-- @m_div{m-button m-primary}
<a href="https://pressio.github.io/pressio-tutorials/html/md_pages_rom_galerkin_default_implicit_custom_types.html">
@m_div{m-big} Implicit Galerkin sample code 2: @m_enddiv
@m_div{m-small} (for arbitrary data types) @m_enddiv
</a> @m_enddiv
 -->