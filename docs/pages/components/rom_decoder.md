
# rom: Arbitrary Decoder

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

A valid decoder is any C++ object whose type meets the following API:

@code{.cpp}
struct Decoder
{
  // these nested typedefs are mandatory because pressio detects them
  using jacobian_type  = /* your type */;
  using fom_state_type = /* your type */;

public:
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
	//  - if romOperand is rank-1, it supports the (i) operator to reference the i-th element
	//  - if romOperand is rank-2, it supports the (i,j) operator to reference the i,j-th element
    // ...
  }

  // if applicable, update the Jacobian for a given state
  template <typename OperandType>
  void updateJacobian(const OperandType & romOperand);

  // return a const reference to the Jacobian matrix
  const jacobian_type & jacobianCRef() const;
};
@endcode

### Requirements

- `fom_state_type`: must be copy constructible

- `jacobian_matrix_type`: must be copy constructible


### Special case: linear decoder

The special case of a linear decoder is natively supported in pressio,
see [this page](./md_pages_components_rom_decoder_linear.html).
