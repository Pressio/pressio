.. role:: raw-html-m2r(raw)
   :format: html


Decoder Abstraction in Pressio
==============================

\todo: complete

A key assumption of projection-based ROMs is to approximate the full-order 
model (FOM) state, @f$y\ *{fom}@f$, as: 
@f[ 
y*\ {fom} = g(y_{rom}) 
@f] 

where @f$y_{rom}@f$ is the reduced state (or generalized coordinates), 
and @f$g@f$ is the mapping between the two. 

Linear
------

If @f$g@f$ is linear, then we can write: 
@f[ 
y\ *{fom} = \phi y*\ {rom} 
@f] 
where @f$\phi@f$ is a matrix (for the time being, assume it constant). 
The Jacobian of the mapping is: 
@f[ 
\frac{d y\ *{fom}}{d y*\ {rom}} = \phi. 
@f] 

Graphically, this corresponds to: 
@image html decoder.png width=20% 

In the C++ library, the class representing the linear decoder is:

@code{.cpp}
pressio::rom::LinearDecoder<...>;
@endcode

Generic decoder API
-------------------

If you want to use your own decoder, you can write a custom one that should meet 
the following API.

@code{.cpp}
struct CustomDecoder
{
  // natives types
  using native_fom_state_type = /\ **/; 
  using native_dense_mat_type = /**\ /; 

  // these nested typedefs are mandatory because pressio detects them
  using jacobian_type  = pressio::containers::DenseMatrix\ :raw-html-m2r:`<native_dense_mat_type>`\ ;
  using fom_state_type = pressio::containers::Vector\ :raw-html-m2r:`<native_fom_state_type>`\ ;

public:
  CustomDecoder(/\ *constuct as needed*\ /){}

  template <typename operand_type>
  void applyMapping(const operand_type & romOperand,\ :raw-html-m2r:`<br>`
                    fom_state_type & fomState) const
  {
    // romOperand: object exposing the (i) operator to reference the i-th element
    // result: a pressio wrapper, so use the data() to get a pointer to native object
    // auto & fomStateNative = *fomState.data();
    // ...
  }

  template <typename rom_state_type>
  void updateJacobian(const rom_state_type & romStateIn)
  {
    // update the jacobian as needed
    // ...
  }

  const jacobian_type & jacobianCRef() const{ return m_jacobian; }

private:
  jacobian_type m_jacobian;
};
@endcode
