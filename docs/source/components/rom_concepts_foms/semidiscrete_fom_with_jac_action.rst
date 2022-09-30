.. include:: ../../mydefs.rst

``SemiDiscreteFomWithJacobianAction``
=====================================

Header: ``<pressio/rom_concepts.hpp>``

.. literalinclude:: ../../../../include/pressio/rom/concepts/fom_semi_discrete_with_jac_action.hpp
   :language: cpp
   :lines: 58-100


The concept ``SemiDiscreteFomWithJacobianAction`` refines
the ``SemiDiscreteFom`` by adding support for the Jacobian evaluation.

Semantic requirements
---------------------

Given an instance ``A`` of type ``T`` and an object ``operand``
of type ``JacobianActionOperandType``,
``SemiDiscreteFomWithJacobianAction<T, JacobianActionOperandType>``
is modeled if it is satisfied, all subsumed concepts are modeled and:

- all methods are blocking, meaning that all temporary
  allocations and operations needed to execute those methods
  are completed and no outstanding work remains upon return

- methods may modify only the non-constant operands.
  Operands that are constant must not be modified.

- ``auto result = A.createApplyJacobianResult(operand)``

  - returns an object with all its "elements" zero initialized


- doing:

  .. code-block:: cpp

     auto ja1 = A.createApplyJacobianResult(operand);
     auto ja2 = A.createApplyJacobianResult(operand);
     // ...
     auto jaN = A.createApplyJacobianResult(operand);

  implies that ``ja1, ja2, ..., jaN`` must be distinct objects,
  and such that any modification to ``ja1`` does not affect any of the others and vice versa.
  In other words, calling ``A.createApplyJacobianResult(operand)`` yields independent instances.

- ``A.applyJacobian(state, operand, evalTime, result)``

  - overwrites ``result`` with the result

  - is equality preserving, i.e. given equal
    inputs ``state, evalTime``, the result remains the same

  - let ``J`` represent the Jacobian which we compute the action of,
    then ``J`` must be the jacobian of the right hand side evaluated for the
    same ``state`` and ``evalTime``. In other words, the Jacobian used for
    computing its action must be mathematically "consistent" with the right hand side.


Syntax-only example
-------------------

.. code-block:: cpp

   class SampleClass
   {
     public:
       using time_type            = double;
       using state_type           = Tpetra::Vector<>; // uses default template parameters
       using right_hand_side_type = state_type;

       right_hand_side_type createRightHandSide() const;
       void rightHandSide(const state_type &     /*state*/,
                          time_type              /*evalTime*/,
			  right_hand_side_type & /*result*/) const;

       Tpetra::MultiVector<> createApplyJacobianResult(const Tpetra::MultiVector<> & operand);

       void applyJacobianResult(const state_type & /*state*/,
                                const Tpetra::MultiVector<> & /*operand*/,
				time_type /*evalTime*/,
                                const Tpetra::MultiVector<> & /*result*/);
   }


Assuming the default scalar type of Tpetra is ``double``,
the class above satisfies: ``static_assert(pressio::rom::SemiDiscreteFomWithJacobianAction<SampleClass,
Tpetra::MultiVector<>>, "");``
