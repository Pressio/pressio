Required ops
============

This page lists all operations used by the nonlinear solvers.
Note that if you use the nonlinear solvers
using *already supported* types (i.e., types which pressio already
knows how to operate on efficiently), then these operations
are already implemented in pressio. If you use *arbitrary data types*,
then all you have to do is to specialize these ops for those types.

:red:`finish  this`

.. code-block:: cpp

   namespace pressio{ namespace ops{

   std::size_t extent(<operand_type> &, int dim);

   <return_type> clone(const <operand_type> & src);

   void deep_copy(<operand_type> & to, const <operand_type> & from);

   void set_zero(ACustomStateType & object);

   <return_type> norm2(const <operand_type> &);

   void update(<operand_type> & v,  scalar_type a,
               const <operand_type> & v1, scalar_type b);

   <return_type> product( for rank2)

   void product( for rank2)

   void scale(...)

   void fill(...)

   <return_type> dot( ... )

   }}
