
# Tutorial: Linear Decoder with arbitrary data structures

@m_class{m-block m-info}

@par Content
This tutorial shows how to create a linear decoder for
data structures NOT supported in pressio. This is the scenario
where you have an application using some arbitrary data types
for which pressio does not know how to operate on.


## Context
A key assumption of projection-based ROMs
is to approximate the full-order
model (FOM) state, @f$y_{fom}@f$, as:
@f[
y_{fom} = g(y_{rom})
@f]

where @f$y_{rom}@f$ is the reduced state, also called
generalized coordinates, and @f$g@f$ is the mapping between the two.
If @f$g@f$ is linear, then we can write:
@f[
y_{fom} = \phi y_{rom}
@f]
where @f$\phi@f$ is a matrix (for the time being, assume it constant).
A linear decoder in pressio implements this linear mapping.
Since it is linear, the Jacobian of the mapping is:
@f[
\frac{d y_{fom}}{d y_{rom}} = \phi.
@f]

## Code
Here we demonstate how to create a linear decoder object for a type that
is NOT know to pressio: this means pressio does not know how to compute
operations on this type, so the user is responsible to pass the ops.
The full tutorial can be found [here](https://github.com/Pressio/pressio-tutorials/blob/master/tutorials/tutorial5.cc).

```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/tutorial5.cc
75:154
```

## The ops struct
In order for pressio to handle the linear mapping, it needs to know
how to operate on @f$\phi@f$. To this end, in the code above,
you need to pass to the `LinearDecoder` constructor an object that
handle that computation.
```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/tutorial5.cc
51:73
```
