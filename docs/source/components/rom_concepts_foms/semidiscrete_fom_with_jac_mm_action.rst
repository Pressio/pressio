.. include:: ../../mydefs.rst

``SemiDiscreteFomWithJacobianAndMassMatrixAction``
==================================================

Header: ``<pressio/rom_concepts.hpp>``

.. literalinclude:: ../../../../include/pressio/rom/concepts/fom_semi_discrete_with_jac_and_mm_action.hpp
   :language: cpp
   :lines: 54-63


The concept ``SemiDiscreteFomWithJacobianAndMassMatrixAction`` refines
the ``SemiDiscreteFomWithJacobianAction`` and ``SemiDiscreteFomWithMassMatrixAction``.

Semantic requirements
---------------------

TBD


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

       Tpetra::MultiVector<> createResultOfJacobianActionOn(const Tpetra::MultiVector<> & operand) const;

       Tpetra::MultiVector<> createResultOfMassMatrixActionOn(const Tpetra::MultiVector<> & operand) const;

       void applyJacobian(const state_type & /*state*/,
                          const Tpetra::MultiVector<> & /*operand*/,
			  time_type /*evalTime*/,
                          const Tpetra::MultiVector<> & /*result*/) const;

       void applyMassMatrix(const state_type & /*state*/,
		            const Tpetra::MultiVector<> & /*operand*/,
			    time_type /*evalTime*/,
			    const Tpetra::MultiVector<> & /*result*/) const;
   }
