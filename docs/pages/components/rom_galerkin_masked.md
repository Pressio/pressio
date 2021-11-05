
# rom: Galerkin: masked problem


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
  class FomReferenceStateType,
  class ProjectorType,
  class MaskerType
  >																				 (1)
ReturnType create_masked_explicit_problem(pressio::ode::StepScheme scheme,
										  const FomSystemType & fomSystem,
										  DecoderType & decoder,
										  const RomStateType & romState,
										  const FomReferenceStateType & fomRefState,
										  const ProjectorType & projector,
										  const MaskerType & masker);

template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class ProjectorType,
  class MaskerType
  >																				 (2)
ReturnType create_masked_implicit_problem(pressio::ode::StepScheme scheme,
										  const FomSystemType & fomSystem,
										  DecoderType & decoder,
										  const RomStateType & romState,
										  const FomReferenceStateType & fomRefState,
										  const ProjectorType & projector,
										  const MaskerType & masker);
```

This function returns an instance of the desired Galerkin problem.

### Parameters and Requirements

- `scheme`:
  - enum value to specify the stepper scheme
  - (1) explicit Galerkin, see [valid enum scheme choices](md_pages_components_ode_steppers_explicit.html)
  - (2) implicit Galerkin, seee [valid enum scheme choices](md_pages_components_ode_steppers_implicit.html)

- `fomSystem`:
  - instance of your FOM adapter specifying the FOM problem <br/>
  - must satisfy one of the APIs suitable for Galerkin, see [API list](./md_pages_components_rom_fom_apis.html)

- `decoder`:
  - decoder object
  - must satify the requirements listed [here](md_pages_components_rom_decoder.html)

- `romState`:
  - ROM state
  - currently, it must be either an Eigen vector or a Kokkos 1D view

- `fomRefState`:
  - your FOM reference state that is used when reconstructing the FOM state
  - must be copy-constructible and the following must be true:<br/>
  ```cpp
  std::is_same<FomReferenceStateType, typename DecoderType::fom_state_type>::value == true
  ```

- `projector`:
  - operator is responsible for projectng the FOM operators onto the reduced space
  - must be a functor with a specific API, see [this page](md_pages_components_rom_galerkin_hypred.html)

- `masker`:
  - functor responsible of "masking" the FOM operators
  - must be a functor with a specific API, see details below


## Masker

\todo: explain

Suppose that `fom_velocity_type` is the type of your FOM velocity,
and that `decoder_jacobian_t` is the type
you use to represent the Jacobian of your decoder function.
To be concrete, for the sake of explanation, `fom_velocity_type` can be for
example a Trilinos
vector type, or a PETSc vector, or any other thing you want.
Similarly for the `decoder_jacobian_t`.

The masker must be a functor as follows:

```cpp
struct ValidMasker
{
  fom_velocity_type createApplyMaskResult(const fom_velocity_type & unmasked_object);

  template<class TimeType>
  void operator()(const fom_velocity_type & unmasked_object,
                  const TimeType time,
				  fom_velocity_type & result);

  // the two methods below are ONLY needed if you are doing implicit time
  decoder_jacobian_type createApplyMaskResult(const decoder_jacobian_type & unmasked_object);

  template<class TimeType>
  void operator()(const decoder_jacobian_type & unmasked_object,
                  const TimeType time,
				  decoder_jacobian_type & result);
};
```
