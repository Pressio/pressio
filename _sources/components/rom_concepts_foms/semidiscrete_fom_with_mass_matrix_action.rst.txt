.. include:: ../../mydefs.rst

``SemiDiscreteFomWithMassMatrixAction``
=======================================

Header: ``<pressio/rom_concepts.hpp>``

.. literalinclude:: ../../../../include/pressio/rom/concepts/fom_semi_discrete_with_mm_action.hpp
   :language: cpp
   :lines: 54-94


The concept ``SemiDiscreteFomWithMassMatrixAction`` refines
the ``SemiDiscreteFom`` by adding support for the mass matrix evaluation.

Semantic requirements
---------------------

Given an instance ``A`` of type ``T`` and an instance ``operand`` of
type ``MassMatrixActionOperandType``,
``SemiDiscreteFomWithMassMatrixAction<T, MassMatrixActionOperandType>``
is modeled if it is satisfied and all of the following are true:

- methods are blocking: all temporary allocations and operations
  needed to execute those methods
  are completed and no outstanding work remains upon return

- methods only modify non-constant arguments, while const arguments are not modified

- ``auto result = A.createResultOfMassMatrixActionOn(operand)`` returns an object
  with all its "elements" zero initialized

- non-aliasing instantation: this means that doing ``auto ja1 = A.createResultOfMassMatrixActionOn(operand);
  auto ja2 = A.createResultOfMassMatrixActionOn(operand);`` implies that ``ja1, ja2`` are distinct objects,
  and such that any modification to ``ja1`` does not affect ``ja2`` and vice versa.
  In other words, calling ``A.createResultOfMassMatrixActionOn(operand)``
  yields **independent, non-aliasing instances**.

- ``A.applyMassMatrix(state, operand, evalTime, result)``

  - overwrites ``result`` with the result of left-applying the mass matrix to ``operand``

  - is equality preserving, i.e. given equal inputs ``state, evalTime``, the result remains the same

  - let ``M`` represent the mass matrix which we compute the action of; then ``M`` must be
    the correct mass matrix for the target problem evaluated at the given
    ``state`` and ``evalTime``. Obivously, if your problem has a constant mass matrix,
    you don't have to recomputed the mass matrix every time.


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

       Tpetra::MultiVector<> createResultOfMassMatrixActionOn(const Tpetra::MultiVector<> & operand) const;

       void applyMassMatrix(const state_type & /*state*/,
		            const Tpetra::MultiVector<> & /*operand*/,
			    time_type /*evalTime*/,
			    const Tpetra::MultiVector<> & /*result*/) const;
   }


..
  Assuming the default scalar type of Tpetra is ``double``,
  the class above satisfies: ``static_assert(pressio::rom::SemiDiscreteFomWithMassMatrixAction<SampleClass,
  Tpetra::MultiVector<>>, "");``
