
# Adapter/system class API for Galerkin ROM

## Continuous-time API

The continuous-time API operates such that the user is responsible
to compute the continuous-time operators, e.g., the velocity,
and pressio assembles the discrete-time operators.


## Basic Galerkin ROM

@m_class{m-code-figure} @parblock
@code{.cpp}
class AdapterGalerkin{
public:
  using scalar_type =
  using state_type =
  using velocity_type =

public:
  velocity_type createVelocity() const;
  void velocity(state, time, velo) const;
};
@endcode
@endparblock


## Discrete-time API
