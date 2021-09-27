
# rom: LSPG: steady default problem


@m_class{m-note m-default}

@parblock
Defined in: `<pressio/rom_lspg.hpp>`

Public namespace: `pressio::rom::lspg`
@endparblock


## Overview


## 1. Creating a problem instance


```cpp
template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType
  >
auto create_default_steady_problem(const FomSystemType &,
								   DecoderType &,
								   const RomStateType &,
								   const FomReferenceStateType &);

template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class PreconditionerType
  >
auto create_default_steady_problem(const FomSystemType &,
								   DecoderType &,
								   const RomStateType &,
								   const FomReferenceStateType &,
								   const PreconditionerType &);
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


## 2.,3. Extract and Solve

```cpp
int main()
{
// we assume the rom_state, decoder, fom_system, fom_reference_state already exist

namespace plspg = pressio::rom::lspg;

auto problem = plspg::create_default_steady_problem(fom_system, decoder,
													rom_state, fom_reference_state);
auto & system = problem.system();

// create nonlinear least-squares solver, for example:
auto nonlinsolver = pressio::ode::create_gauss_newton(system, rom_state, ...);
nonlinsolver.solve(system, rom_state);
//...
}
```
