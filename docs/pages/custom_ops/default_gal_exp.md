
# Custom Ops for Default Explicit-time Galerkin

This page explains which operations you need to provide
to pressio to construct and run Galerkin with explicit time stepping
when you are using FOM types that are *not* natively supported in pressio.

We remark that when you types already known to pressio,
you don't need to provide any custom ops.

```cpp
template <typename scalar_t>
struct CustomOps // or whatever name you want
{
  void deep_copy(TypeOfYourForFomState & dest,
				 const TypeOfYourForFomState & src) const
  {
    // deep copy content of src into dest
  }

  void set_zero(TypeOfYourForFomState & vec) const
  {
	 // self explanatory
  }

  void axpy(scalar_t alpha,
		    const TypeOfYourForFomState & x,
			TypeOfYourForFomState & y) const
  {
    // compute y = y + alfa * x
  }

  /*
  compute: y = beta * y + alpha*A*x
  - y is your application vector type
  - A is the basis matrix
  - x is the rom state
  */
  template <typename x_t>
  void product(::pressio::nontranspose mode,
			   const scalar_t alpha,
	           const TypeUsedForBasisMatrix & A,
			   const x_t & x,
			   const scalar_t beta,
			   TypeOfYourForFomState & y) const
  {
    // when running on CPU, here you can assume that
	// x is subscriptable on host as: x(i)
  }

  /*
  compute: y = beta * y + alpha*A^T*x
  - y is the rom state
  - A is the basis matrix
  - x is your application vector type
  */
  template <typename y_t>
  void product(::pressio::transpose mode,
			   const scalar_t alpha,
			   const TypeUsedForBasisMatrix & A,
			   const TypeOfYourForFomState & x,
			   const scalar_t beta,
			   y_t & y) const
  {
    // when running on CPU, here you can assume that
	// y is subscriptable on host as: y(i)
  }
};
```

## Where does one need to use this?

An object of the class template above needs to be
passed to the constructor of the Galerkin problem as
explained in [this tutorial](./md_pages_tutorials_tutorial2.html).
