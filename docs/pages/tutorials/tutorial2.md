
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
int main(int argc, char *argv[])
{
  // *** define some types ***
  // here we assume your FOM application uses an Eigen vector for the state
  // and an Eigen matrix as the type for the Jacobian
  using scalar_t	       = double;
  using native_fom_state_t = std::vector<scalar_t>;
  using native_phi_t	   = std::vector<std::vector<scalar_t>>;

  // the wrapped types
  // what happens in pressio: std::vector is treated as unknwon type by pressio
  // so effectively pressio::containers::Vector is labeled as an "arbitrary" type
  using fom_state_t	   = pressio::containers::Vector<native_fom_state_t>;
  using decoder_jac_t  = pressio::containers::DenseMatrix<native_phi_t>;

  // *** fill phi ***
  // create a native phi with and fill with ones
  native_phi_t phiNative(6);
  for (auto & iRow : phiNative){
    iRow.resize(2, 1.);
  }

  // *** construct decoder  ***
  using ops_t = MyOps<scalar_t>;
  using decoder_t = pressio::rom::LinearDecoder<decoder_jac_t, fom_state_t, ops_t>;
  // Need to pass the native phi (here we assume the native type is copy-constructible)
  // and an object that knows how to compute the operations (see next section)
  ops_t ops;
  decoder_t decoder(phiNative, ops);

  // *** construct reduced state  ***
  // typically, pressio reduced states for ROMs use Eigen or Kokkos (if enabled)
  using rom_state_t = pressio::containers::Vector<Eigen::VectorXd>;
  rom_state_t yRom(2);
  // set yRom = 2.
  pressio::ops::fill(yRom, 2.);

  // *** apply mapping ***
  fom_state_t yFom(6);
  decoder.applyMapping(yRom, yFom);

  // *** check solution ***
  // yFom should be = [4. 4. .... 4.]
  for (auto i=0; i<6; ++i)
    std::cout << "i= " << i
	      << ", yFom(i) = "  << yFom(i)
	      << ", expected = " << 4.
	      << "\n";

  return 0;
}
```

## The ops struct
In order for pressio to handle the linear mapping, it needs to know
how to operate on @f$\phi@f$. To this end, in the code above,
you need to pass to the `LinearDecoder` constructor an object that
handle that computation.
```cpp

template <typename scalar_t>
struct MyOps
{
  // z = beta*z + alpha * A * x
  // where x is something that is subscritable as x(i)
  template< typename x_t>
  void product(pressio::nontranspose,
	       scalar_t alpha,
	       const std::vector<std::vector<scalar_t>> & A,
	       const x_t & x,
	       scalar_t beta,
	       std::vector<scalar_t> & z) const
  {
    // obviously not efficient, just for demonstration
    for (std::size_t i=0; i<A.size(); ++i)
    {
      z[i] += beta*z[i];
      for (std::size_t j=0; j<A[i].size(); ++j){
	   z[i] += alpha*A[i][j]*x(j);
      }
    }
  }
};
```
