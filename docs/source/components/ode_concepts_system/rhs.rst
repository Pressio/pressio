.. include:: ../../mydefs.rst

``SystemWithRhs``
=================

.. literalinclude:: ../../../../include/pressio/ode/concepts/system_rhs.hpp
   :language: cpp
   :lines: 52-93


Semantic requirements
---------------------

:red:`finish`


Syntax-only example
-------------------

.. code-block:: cpp

   class SampleClass
   {
     public:
       using independent_variable_type = double;
       using state_type                = Eigen::VectorXd;
       using right_hand_side_type      = Eigen::VectorXd;

       state_type createState() const;
       right_hand_side_type createRightHandSide() const;
       void operator()(const state_type &        /*state*/,
                       independent_variable_type /*ivEval*/,
		       right_hand_side_type &    /*rhsToOverwrite*/) const;
   }



..
   The concept is modeled only if all of the following hold:

   - *determined*: the "dimension" of the right hand side must be equal to the
     dimension of the state, i.e. # of equations is same as # of unknowns

   - *non aliasing operators*: given the following:

     .. code-block:: cpp

	auto s1 = A.createState();
	auto s2 = A.createState();
	auto r1 = A.createRightHandSide();
	auto r2 = A.createRightHandSide();

     ``r1``, ``r2`` must be distinct objects, ``std::addressof(r1) != std::addressof(r2)``,
     and such that any modification to ``r1`` does not affect ``r2``, and similarly
     for ``s1`` and ``s2``

   - *blocking operations*: all methods are blocking, meaning that all temporary
     allocations and operations are complete before the methods return and not outstanding work remains

   - *equality preserving*: given ``A`` an object of type `T`, calling ``A.rightHandSide(...)``
     with equal inputs yields equal outputs

   - *const correctness*: methods may modify only the non-constant operands.
     Operands that are constant must not be modified.
