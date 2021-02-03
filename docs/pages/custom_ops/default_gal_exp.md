
# Custom Ops for Default Explicit-time Galerkin

```cpp
template <typename scalar_t>
struct CustomOps
{

  //  y = beta * y + alpha*A*x
  template <typename x_t>
  void product(::pressio::nontranspose mode,
			   const scalar_t alpha,
	           const dec_jac_t & A,
			   const x_t & x,
			   const scalar_t beta,
			   pressio::apps::arbds::Vector<scalar_t> & y) const
  {
    // x is subscriptable like a regular array, e.g. you can do x[i] or x(i)
  }

  //  y = beta * y + alpha*A^T*x
  template <typename y_t>
  void product(::pressio::transpose mode,
	       const scalar_t alpha,
	       const dec_jac_t & A,
	       const pressio::apps::arbds::Vector<scalar_t> & x,
	       const scalar_t beta,
	       y_t & y) const
  {
    // y is subscriptable like a regular array
  }

  void deep_copy(pressio::apps::arbds::Vector<scalar_t> & to,
  		 const pressio::apps::arbds::Vector<scalar_t> & from) const
  {
    // here you need do deep copy from -> to
  }

  void set_zero(pressio::apps::arbds::Vector<scalar_t> & vec) const
  {
	 // self explanatory
  }

  void axpy(scalar_t alpha,
	    const pressio::apps::arbds::Vector<scalar_t> & x,
	    pressio::apps::arbds::Vector<scalar_t> & y) const
  {
    // compute y = y + alfa * x
  }
};
};
```
