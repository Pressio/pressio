
# rom: Galerkin: default problem


Defined in: `<pressio/rom_galerkin.hpp>`

Public namespace: `pressio::rom`, `pressio::rom::galerkin`.


## Overview
bla bla


## API

```cpp
template<
  class StepperTag,
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType
  >
auto create_default_problem(const FomSystemType & fomSysObj,
                            DecoderType & decoder,
                            const RomStateType & stateIn,
                            const FomReferenceStateType & fomRef);
```

where `StepperTag` is a tag type from the ode to specify which time scheme to use.
This function returns an instance of the desired Galerkin problem.

## Parameters and Requirements

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
  - currently, it must be either an Eigen or Kokkos type

- `FomReferenceStateType`:
  - your FOM reference state that is used when reconstructing the FOM state
  - must be copy-constructible and the following must be true:<br/>
  ```cpp
  std::is_same<FomReferenceStateType, typename DecoderType::fom_state_type>::value == true
  ```

## What is a problem and how to use it?

When you create a problem object, pressio behind
the scenes checks types, and instantiates everything that is needed
to define a problem. Practically speaking, at the lowest-level,
a Galerkin problem reduces to creating a "custom" stepper to advance in time.
This is exactly how pressio implements this: when you create a Galerkin
problem, pressio behind the scenes creates this custom stepper
object that you can use. You don't know to know or rely on
the details of the stepper or how this is done,
because these are problem-dependent and we reserve the right
to change this in the future.

All you need to know to use a problem is that it has the following interface:

```cpp
class DefaultGalerkinProblem
{
public:
  // returns the underlying stepper to use to solve the problem
  auto & stepper();

  // const ref to the object knowing how to reconstruct a FOM state
  const auto & fomStateReconstructor() const;
};
```

where the stepper is an [explicit stepper](md_pages_components_ode_steppers_explicit.html)
if you used an explicit `StepperTag` or an [implicit stepper](md_pages_components_ode_steppers_implicit.html)
if you used an implicit `StepperTag`.
Once you have the stepper, you can then use it to solve your problem like
you would with any other stepper object.
An example usage is as follows:

```cpp
int main()
{
// we assume the rom_state, decoder, fom_system, fom_reference_state already exist

namespace pode = pressio::ode;
namespace pgal = pressio::rom::galerkin;

using ode_tag  = pode::ForwardEuler;
auto problem   = pgal::create_default_problem<ode_tag>(fom_system, decoder, rom_state, fom_reference_state);
auto & stepper = problem.stepper();

pressio::ode::advance_n_steps_and_observe(stepper, rom_state, /* any other args */);
}
```
