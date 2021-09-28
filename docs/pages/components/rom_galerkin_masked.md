
# rom: Galerkin: masked problem


@m_class{m-note m-default}

@parblock
Defined in: `<pressio/rom_galerkin.hpp>`

Public namespace: `pressio::rom::galerkin`
@endparblock


## Overview

Recall from [this page](md_pages_components_rom_galerkin.html),
that using a pressio Galerkin problem involves three steps:

1. *create*: you create an instance of a "Galerkin problem"

2. *extract*: you extract the underlying stepper object owned by the problem

3. *solve*: you use the stepper to solve in time the Galerkin problem


<br/>
___
<br/>


## 1. Creating a problem instance

```cpp
template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class ProjectorType,
  class MaskerType
  >
ReturnType create_masked_explicit_problem(pressio::ode::StepScheme,
										  const FomSystemType &,
										  DecoderType &,
										  const RomStateType &,
										  const FomReferenceStateType &,
										  const ProjectorType &,
										  const MaskerType &);

template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class ProjectorType,
  class MaskerType
  >
ReturnType create_masked_implicit_problem(pressio::ode::StepScheme,
										  const FomSystemType &,
										  DecoderType &,
										  const RomStateType &,
										  const FomReferenceStateType &,
										  const ProjectorType &,
										  const MaskerType &);
```

This function returns an instance of the desired Galerkin problem.

### Parameters and Requirements

- `StepScheme`:
  - must be one of the explicit or implicit enum values supported in `pressio::ode`

- `FomSystemType`:
  - your adapter class type specifying the FOM problem
  - must satisfy one of the APIs suitable for Galerkin, see [API list](./md_pages_components_rom_fom_apis.html)

- `DecoderType`:
  - decoder class type
  - must satify the requirements listed [here](md_pages_components_rom_decoder.html)

- `RomStateType`:
  - ROM state type
  - currently, it must be either an Eigen vector or a Kokkos 1D view

- `FomReferenceStateType`:
  - your FOM reference state that is used when reconstructing the FOM state
  - must be copy-constructible and the following must be true:<br/>
  ```cpp
  std::is_same<FomReferenceStateType, typename DecoderType::fom_state_type>::value == true
  ```

- `ProjectorType`:
  - an operator is responsible for projectng the FOM operators onto the reduced space
  - must be a functor with a specific API, see [this page](md_pages_components_rom_galerkin_hypred.html)

- `MaskedType`:
  - an functor responsible of "masking" the FOM operators
  - must be a functor with a specific API, see details below


### Masker

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
