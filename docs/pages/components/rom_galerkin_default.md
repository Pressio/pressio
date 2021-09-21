
# rom: Galerkin: default problem


Defined in: `<pressio/rom_galerkin.hpp>`

Public namespace: `pressio::rom::galerkin`


## Overview

At a high level, using a Galerkin problem involces three steps:
1. *create*: you create an instance of a "default Galerkin problem"
2. *extract*: you extract the underlying stepper object owned by the problem
3. *solve*: you use the stepper to solve in time the Galerkin problem


You should now pause and think for a second about the steps above.
What does a stepper have to do with a Galerkin ROM?
The answer is that practically speaking, at the lowest-level,
a Galerkin problem can be reduced to simply a "custom" stepper to advance in time.
This is exactly how pressio implements this and the reason why a Galerkin
problem contains a stepper object inside: when you create the
problem, pressio creates the appropriate custom stepper
object that you can use. You don't need to know how this is done,
or rely on the details, because these are problem- and implementation-dependent,
and we reserve the right to change this in the future.


## 1. Creating a problem instance

```cpp
template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType
  >
auto create_default_explicit_problem(pressio::ode::StepScheme,
								     const FomSystemType &,
									 DecoderType &,
									 const RomStateType &,
									 const FomReferenceStateType &)

template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType
  >
auto create_default_implicit_problem(pressio::ode::StepScheme,
								     const FomSystemType &,
									 DecoderType &,
									 const RomStateType & ,
									 const FomReferenceStateType &)
```

This function returns an instance of the desired Galerkin problem.

### Parameters and Requirements

- `StepScheme`:
  - enum value to specify the stepper scheme, see [explicit choices](md_pages_components_ode_steppers_explicit.html) and [implicit choices](md_pages_components_ode_steppers_implicit.html)

- `FomSystemType`:
  - your adapter class type specifying the FOM problem. <br/>
  - Must satisfy one of the APIs suitable for Galerkin, see [API list](./md_pages_components_rom_fom_apis.html)

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

### Galerkin Problem class API

A problem meets the following interface:

```cpp
class DefaultGalerkinProblem
{
public:
  using traits = /* nested typedef with trait class */;

  // returns the underlying stepper to use to solve the problem
  auto & stepper();

  // const ref to the object knowing how to reconstruct a FOM state
  const auto & fomStateReconstructor() const;
};
```

The `stepper` method is, practically, what you would use
to retrieve the stepper and then use it to solve the problem.
The stepper method returns a non-const reference to an
[explicit stepper](md_pages_components_ode_steppers_explicit.html)
if, when you created the problem, you used an explicit `StepperTag`,
or an [implicit stepper](md_pages_components_ode_steppers_implicit.html)
if you use an implicit `StepperTag`.
Once you retrieve the stepper, you can then use it like
you would with any other stepper object (more on this below).

As almost every important type in pressio, a Galerkin problem
too has traits:

```cpp
auto problem = create_default_problem(...);
using traits = pressio::Traits<decltype(problem)>;

// for the explicit case, one can access the following traits:
typename traits::fom_system_type;
typename traits::scalar_type;
typename traits::decoder_type;
typename traits::decoder_jac_type;
typename traits::galerkin_state_type;
typename traits::galerkin_velocity_type;
typename traits::stepper_type;

// for the implicit case one has:
typename traits::fom_system_type;
typename traits::scalar_type;
typename traits::decoder_type;
typename traits::decoder_jac_type;
typename traits::galerkin_state_type;
typename traits::galerkin_residual_type;
typename traits::galerkin_jacobian_type;
typename traits::stepper_type;
```

## 2.,3. Extract and Solve

### Explicit Case
An example usage for explicit stepper is as follows:

```cpp
int main()
{
// we assume the rom_state, decoder, fom_system, fom_reference_state already exist

namespace pode = pressio::ode;
namespace pgal = pressio::rom::galerkin;

const auto scheme = pdoe::StepScheme:ForwardEuler;
auto problem = pgal::create_default_explicit_problem(scheme, fom_system, decoder,
													 rom_state, fom_reference_state);
auto & stepper = problem.stepper();

pressio::ode::advance_n_steps_and_observe(stepper, rom_state, /* any other args */);
}
```

### Implicit Case
An example usage for implicit stepper is as follows:

```cpp
int main()
{
// we assume the rom_state, decoder, fom_system, fom_reference_state already exist

namespace pode = pressio::ode;
namespace pgal = pressio::rom::galerkin;

const auto scheme = pdoe::StepScheme:BDF1;
auto problem = pgal::create_default_implicit_problem(scheme, fom_system, decoder,
													 rom_state, fom_reference_state);
auto & stepper = problem.stepper();

!!!!!!!! fill !!!!!!!!!

pressio::ode::advance_n_steps_and_observe(stepper, rom_state, /* any other args */);
}
```
