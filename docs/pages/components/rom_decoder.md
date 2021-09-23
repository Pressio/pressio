
# rom: Decoder

## Overview

A key assumption of projection-based ROMs is to approximate
the full-order model (FOM) state, @f$y_{fom}@f$, as:
@f[
y_{fom} = g(y_{rom})
@f]

where @f$y_{rom}@f$ is the reduced state (or generalized coordinates),
and @f$g@f$ is the mapping between the two.

Note that there is no explicit constraint on what the mapping is, it can be anything.

## Decoder Concept

In pressio, a valid decoder is any C++ object whose type meets the following API:

@code{.cpp}
class Decoder
{
public:
  // these nested typedefs are mandatory because pressio detects them
  using jacobian_type  = /* your type */;
  using fom_state_type = /* your type */;

  template <class OperandType>
  void applyMapping(const OperandType & romOperand,
                    fom_state_type & fomState) const
  {
    // romOperand:
	//  typically, this is the ROM state (or generalized coordinates),
	//  but this is not always true: in some cases, e.g., WLS, it can be an expression.
	//  In general, we advise to keep it as a template.
	//
	//  To use romOperand, you need to know that:
	//  - rank-1 romOperand: supports the (i) operator to reference the i-th element
	//  - rank-2 romOperand: supports the (i,j) operator to reference the i,j-th element
    // ...
  }

  // if applicable, update the Jacobian for a given state
  template <typename OperandType>
  void updateJacobian(const OperandType & romOperand);

  // return a const reference to the Jacobian object
  const jacobian_type & jacobianCRef() const;
};
@endcode

### Requirements

- `fom_state_type`: must be copy constructible

- `jacobian_type`: must be copy constructible

<br/>

## Special case: linear decoder


@m_class{m-frame}

@parblock
Defined in header: `<pressio/rom_decoder.hpp>`

Public namespace: `pressio::rom`
@endparblock


A linear decoder is a mapping of the form:
@f[
y_{fom} = \phi y_{rom}
@f]

where @f$\phi@f$ is the Jacobian matrix (for the time being, assume it constant). <br/>
Pressio offers a class for this abstraction:

@code{.cpp}
template<class FomStateType, class JacobianType>
auto create_time_invariant_linear_decoder(JacobianType &&);
@endcode

where:
- `FomStateType`: type of your FOM state, must be copy constructible

- `JacobianType`: type of the decoder's Jacobian, must be copy constructible and move constructible

Obviously, the returned linear decoder object meets the decoder concept discussed above.

### Example usage for supported types

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

### Example usage for custom types

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
