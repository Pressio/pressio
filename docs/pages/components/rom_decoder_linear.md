
# rom: Linear Decoder

Defined in header: `<pressio/rom_decoder.hpp>`

Public namespace: `pressio::rom`

## Overview

A linear decoder is a mapping of the form:
@f[
y_{fom} = \phi y_{rom}
@f]

where @f$\phi@f$ is the Jacobian matrix (for the time being, assume it constant). <br/>
Pressio offers a class for this abstraction
as described below.

@code{.cpp}
template<class FomStateType, class JacobianMatrixType>
auto create_time_invariant_linear_decoder(JacobianMatrixType && jac_matrix);
@endcode

The returned linear decoder object meets
the [decoder concept](./md_pages_components_rom_decoder.html).

### Parameters

- `FomStateType`: data type of your FOM state

- `JacobianMatrixType`: data type your Jacobian matrix


### Requirements

- `FomStateType`: must be copy constructible

- `JacobianMatrixType`: must be copy constructible and move constructible.


## Example usage for supported types

When using data types supported in [pressio ops](./md_pages_components_ops.html), an example usage is as follows:

```cpp
#include "pressio/type_traits.hpp"
#include "pressio/rom_decoder.hpp"
int main()
{
  namespace prom = pressio::rom;

  // assuming that:
  // all proper initialization has been done

  using fom_state_type = Tpetra::Vector<>;
  using matrix_type    = Tpetra::MultiVector<>;

  matrix_type matJ(/* construct as needed */);
  auto decoder = prom::create_time_invariant_linear_decoder<fom_state_type>(matJ);
}
```

## Example usage for custom types

When using custom data types not supported in [pressio ops](./md_pages_components_ops.html),
you need to provide specializations of a trait class and certain operations
and make them "visible" to the compiler to find them and such that pressio can operate on your data.
For the sake of explanation, suppose that you use `double` as value type,
`MyCustomVector<double>` for the FOM state, and `MyCustomMatrix<double>` for the Jacobian matrix.
Then you would need to do something like this:

@code{.cpp}

#include "pressio/type_traits.hpp"

namespace pressio{

template<class ScalarType>
struct Traits<pressiotests::MyCustomVector<ScalarType>>{
  using scalar_type = ScalarType;
};

template<class ScalarType>
struct Traits<pressiotests::MyCustomMatrix<ScalarType>>{
  using scalar_type = ScalarType;
};

namespace ops{

template<class OperandType, class ScalarType>
void product(pressio::nontranspose,
             ScalarType alpha,
             const pressiotests::MyCustomMatrix<ScalarType> & A,
             const OperandType & x,
             ScalarType beta,
             pressiotests::MyCustomVector<ScalarType> & z)
{
  // z = beta*z + alpha * A * x
  // you need to compute a standard gemv, but all you know about x
  // is that its i-th element can be retrieved as x(i)
}

}}//end namespace pressio::ops

#include "pressio/rom_decoder.hpp"
int main()
{
  namespace prom = pressio::rom;

  // assuming that:
  // all proper initialization has been done

  using fom_state_t = MyCustomVector<double>;
  using matrix_t    = MyCustomMatrix<double>;

  matrix_t matJ(/* construct as needed */);
  auto decoder = prom::create_time_invariant_linear_decoder<fom_state_type>(matJ);
}
@endcode
