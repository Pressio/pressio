.. include:: ../../mydefs.rst

``SemiDiscreteFom``
===================

Header: ``<pressio/rom_concepts.hpp>``

.. literalinclude:: ../../../../include/pressio/rom/concepts/fom_semi_discrete.hpp
   :language: cpp
   :lines: 54-96

Semantic requirements
---------------------

Given an instance ``A`` of type ``T``, ``SemiDiscreteFom<T>``
is modeled if it is satisfied and all of the following are true:

- methods are blocking: all temporary allocations and operations
  needed to execute those methods
  are completed and no outstanding work remains upon return

- methods only modify non-constant arguments, while const arguments are not modified

- ``auto rhs = A.createRightHandSide();`` returns an object
  with all its "elements" zero initialized

- non-aliasing instantation: this means that doing ``auto rhs1 = A.createRightHandSide();
  auto rhs2 = A.createRightHandSide();`` implies that ``rhs1, rhs2`` must be distinct objects,
  and such that any modification to ``rhs1`` does not affect ``rhs2`` and viceversa.
  In other words, calling ``A.createRightHandSide()``
  yields **independent, non-aliasing instances**.

- ``A.rightHandSide(state, evalTime, rhs)`` overwrites ``rhs`` with the result,
  and is equality preserving, i.e. given equal inputs ``state, evalTime``, the result
  written to ``rhs`` remains the same


Syntax-only example
-------------------

.. code-block:: cpp

   class SampleClass
   {
     public:
       using time_type            = double;
       using state_type           = Tpetra::Vector<double, /*whatever else*/>;
       using right_hand_side_type = state_type;

       right_hand_side_type createRightHandSide() const;
       void rightHandSide(const state_type &     /*state*/,
                          time_type              /*evalTime*/,
			  right_hand_side_type & /*result*/) const;
   }
