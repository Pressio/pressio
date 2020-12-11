
# Tutorial: Linear Decoder

@m_class{m-block m-info}

@par Content
This tutorial shows how to create a linear decoder for
data structures already supported in pressio.


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
Here we demonstate how to create a linear decoder
object using Eigen types.
The full tutorial can be found [here](https://github.com/Pressio/pressio-tutorials/blob/master/tutorials/tutorial4.cc).

```cpp
int main(int argc, char *argv[])
{
  // *** define some types ***
  // here we assume your FOM application uses an Eigen vector for the state
  // and an Eigen matrix as the type for the Jacobian
  using native_fom_state_t = Eigen::VectorXd;
  using native_phi_t	   = Eigen::MatrixXd;
  // the wrapped types
  using fom_state_t	  = pressio::containers::Vector<native_fom_state_t>;
  using decoder_jac_t = pressio::containers::DenseMatrix<native_phi_t>;

  // *** fill phi ***
  // for simplicity, create a native phi with 10 rows and 3 columns
  // and fill with ones
  native_phi_t phiNative(6, 2);
  phiNative.setConstant(1.);

  // *** construct decoder  ***
  using decoder_t = pressio::rom::LinearDecoder<decoder_jac_t, fom_state_t>;
  // here we demonstrate moving phiNative to avoid a deep copy, but one can
  // also do `decoderObj(phi)` which implies a copy of the matrix.
  // Obviously, if the matrix is large avoiding a copy is useful.
  decoder_t decoder(std::move(phiNative));

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
  for (auto i=0; i<yFom.extent(0); ++i)
	std::cout << "i= " << i
		  << ", yFom(i) = "  << yFom(i)
		  << ", expected = " << 4.
		  << "\n";

  return 0;
}
```

## Comments
For other types (e.g., Trilinos) already known to pressio (which means
pressio knows how to operate on), it would work similarly.
If you work with an arbitrary type (i.e. one for which pressio
does not know how to operate on), see [this tutorial](https://pressio.github.io/pressio/html/md_pages_tutorials_tutorial2.html).
