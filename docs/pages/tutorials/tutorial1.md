
# Tutorial: Linear Decoder

@m_class{m-block m-info}

@par
This tutorial shows how to create a pressio linear decoder object.


## Context
A key assumption of projection-based ROMs is to approximate the full-order
model (FOM) state, @f$y_{fom}@f$, as:
@f[
y_{fom} = g(y_{rom})
@f]

where @f$y_{rom}@f$ is the reduced state (or generalized coordinates),
and @f$g@f$ is the mapping between the two.

If @f$g@f$ is linear, then we can write:
@f[
y_{fom} = \phi y_{rom}
@f]
where @f$\phi@f$ is a matrix (for the time being, assume it constant).
The Jacobian of the mapping is:
@f[
\frac{d y_{fom}}{d y_{rom}} = \phi.
@f]

Graphically, this corresponds to:
@image html tut_f1.png width=35%

*A linear decoder in pressio implements this linear mapping.*


## Types Description

* @f$y_{fom}@f$: this is the FOM state and, therefore, it is a data structure of your application.
For example, if you are using Trilinos/Epetra, the FOM state can be, e.g., an Epetra vector.
Typically, @f$y_{fom}@f$ is a large (distributed) vector.

* @f$\phi@f$: this is the matrix of the linear mapping (e.g. POD modes).
Should be a data structure from your application.
For example, if you are using Trilinos/Epetra, the FOM state can be, e.g., an Epetra MultiVector.
Typically, @f$\phi@f$ is a large (distributed) matrix.

* @f$y_{rom}@f$: this is the ROM state.
ROM data strutures generally involve small, dense operators that fit well on a single node.
In pressio, *regardless of what the FOM types are*, the ROM operators/data structures
are either Eigen (always enabled) or Kokkos (optional dependency) types.



## Scenario A: your FOM types are natively supported in pressio

This case refers to data types which pressio knows how to manipulate and operate on.
Examples include vector and matrix classes in Eigen, Epetra/Tpetra in Trilinos, or Kokkos views.
What do we mean by *natively supported*? If you try to use a pressio
functionality/class using data structure types that are already supported/known
to pressio, pressio *detects* them, and automatically uses the (best)
native kernels to perform computations.
The full list of supported data structures types can be found [here](\todo).

Here, for demonstration purposes, we pretend the FOM uses Eigen types too.
For other FOM types already known to pressio it would work similarly.
The full tutorial can be found [here](https://github.com/Pressio/pressio-tutorials/blob/master/tutorials/tutorial4.cc).

```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/tutorial4.cc
51:104
```


## Scenario B: your FOM types are NOT natively supported in pressio

This case refers to data types which pressio does *not* know how to manipulate and operate on.
What do we mean by *not natively supported*? If you try to use a pressio
functionality/class usign a data structure type that is **NOT** already supported/known
to pressio, pressio *detects/labels* it as an *arbitrary* type, and
you **have to** provide pressio with kernels to operate on these types.

Here, for demonstration, we pretend that:
* @f$y_{fom}@f$ is a `std::vector<>`
* @f$\phi@f$ is a `std::vector<std::vector<>>`

Note that the code below would work for any other "arbitrary'' (not known to pressio) type,
whether they are distributed or not.
The full tutorial can be found [here](https://github.com/Pressio/pressio-tutorials/blob/master/tutorials/tutorial5.cc).

```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/tutorial5.cc
75:137
```

### The ops struct
In order for pressio to handle the linear mapping, it needs to know
how to operate on @f$\phi@f$. To this end, in the code above,
you need to pass to the `LinearDecoder` constructor an object to handle that computation.
To compute the mapping, pressio will call the `product` method.
```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/tutorial5.cc
51:73
```

## Comments

We are constantly working increasing support in pressio for more external libraries.
If you application types are not supported but you would like them to be,
you can file an [issue](https://github.com/Pressio/pressio/issues) to request it.
