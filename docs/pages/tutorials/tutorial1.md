
# Tutorial: Linear Decoder

@m_class{m-block m-info}

@par
This tutorial shows how to create a linear decoder object for
data structures already supported in pressio.


## Context
A key assumption of projection-based ROMs is to approximate the full-order
model (FOM) state, @f$y_{fom}@f$, as:
@f[
y_{fom} = g(y_{rom})
@f]

where @f$y_{rom}@f$ is the reduced state (or generalized coordinates),
and @f$g@f$ is the mapping between the two. <br/>
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
Here we demonstate how to create a linear decoder
object using Eigen types.
The full tutorial can be found [here](https://github.com/Pressio/pressio-tutorials/blob/master/tutorials/tutorial4.cc).

```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/tutorial4.cc
51:122
```

## Comments
For other types (e.g., Trilinos) already known to pressio (which means
pressio knows how to operate on), it would work similarly.
If you work with an arbitrary type (i.e. one for which pressio
does not know how to operate on),
see [this tutorial](tutorial1udops_8md.html).
