
# rom: Galerkin: hyper-reduced problem


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
  class ProjectorType
  >
auto create_hyperreduced_explicit_problem(pressio::ode::StepScheme,
										  const FomSystemType &,
									      DecoderType &,
									      const RomStateType &,
									      const FomReferenceStateType &,
									      const ProjectorType &);

template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class ProjectorType
  >
auto create_hyperreduced_implicit_problem(pressio::ode::StepScheme,
										  const FomSystemType &,
									      DecoderType &,
									      const RomStateType &,
									      const FomReferenceStateType &,
									      const ProjectorType &);
```
This function returns an instance of the desired Galerkin problem.

### Parameters and Requirements

- `StepScheme`:
  - enum value to specify stepper

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
  - must be a functor with a specific API, see details below

### Problem class API

An instance of the hyper-reduced Galerkin problem meets the same API
as the [default problem](md_pages_components_rom_galerkin_default.html).

A hyper-reduced Galerkin problem has these traits:

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
typename traits::projector_type;
typename traits::stepper_type;

// for the implicit case one has:
typename traits::fom_system_type;
typename traits::scalar_type;
typename traits::decoder_type;
typename traits::decoder_jac_type;
typename traits::galerkin_state_type;
typename traits::galerkin_residual_type;
typename traits::galerkin_jacobian_type;
typename traits::projector_type;
typename traits::stepper_type;
```

### Projector

\todo: explain what projectos is for since it is critical.

The projector must be a functor as follows:

```cpp
struct ValidProjector
{
  template<class operand_type, class TimeType, class ResultType>
  void operator()(const operand_type & operand,
                  const TimeType time,
				  ResultType & result) const
  {
    // apply the projection to the operand, store into result
  }
};
```

Note that here the `operator()` is templated on the operand
because the projector can be called with different arguments depending
on whether you are selected an explicit or implicit Galerkin problem:
- for explicit: the projector *always* receives as operand argument
an instance of your FOM velocity
- for implicit: the operand is an instance of your FOM velocity
as well as an instance of the decoder's Jacobian.

You are responsible to handle the various cases.
Obviously, you are not required to template the operand.
If you are working with specific types, you can simply specialize it.
For example, suppose that you are working with an application
that uses for the velocity a type named `FomVelocityType`,
and for the decoder's Jacobian you are using a type named `ADenseMatrixType`.
In such case, a valid projector can be written as:

```cpp
struct ValidProjector
{
  template<class TimeType, class ResultType>
  void operator()(const FomVelocityType & operand,
                  const TimeType time,
				  ResultType & result) const
  {
    // apply the projection to the operand, store into result
  }

  template<class TimeType, class ResultType>
  void operator()(const ADenseMatrixType & operand,
                  const TimeType time,
				  ResultType & result) const
  {
    // apply the projection to the operand, store into result
  }
};
```

The compiler will select the best match.
Note that the result is NOT templated.
In general, you can safely assume that `result` is indexable as `(i,...)`.
\todo describe this more
