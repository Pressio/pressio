
# FOM Adapter API for Galerkin ROM

\todo discuss what is Galerkin ROM, the variants of the problem and what we support.

\todo: discuss that the types exposed must be compatible and make sense.

## Basic Galerkin ROM

@m_class{m-code-figure} @parblock
@code{.cpp}
class Adapter{
public:
  using scalar_type =
  using state_type =
  using velocity_type =
  using dense_matrix_type =

public:
  velocity_type createVelocity() const;
  void velocity(state, time, velo) const;
};
@endcode
@endparblock


## Discrete-time API

The discrete-time API for Galerkin is the same as for all methods,
and you can see it [here](./md_pages_adapter_apis_adapter_discrete_time_api.html).
