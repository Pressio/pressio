
# rom: Galerkin: masked problem

Defined in: `<pressio/rom_galerkin.hpp>`

Public namespace: `pressio::rom`, `pressio::rom::galerkin`.


## 1. Creating a problem instance

```cpp
template<
  class StepperTag,
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class ProjectorType,
  class MaskerType
  >
auto create_masked_problem(const FomSystemType & fomSysObj,
                           DecoderType & decoder,
                           const RomStateType & stateIn,
                           const FomReferenceStateType & fomRef,
						   const ProjectorType & probjector,
						   const MaskerType & masker);
```

where `StepperTag` is a tag type from the ode to specify which time scheme to use.
This function returns an instance of the desired Galerkin problem.

### Parameters and Requirements

- `StepperTag`:
  - tag type to specify the time scheme
  - must be one of the explicit or implicit tag types supported in `pressio::ode`

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
  - an operator responsible of "masking" the FOM operators
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
