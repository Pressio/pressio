.. include:: ../../mydefs.rst

``SteadyFomWithJacobianAction``
===============================

Header: ``<pressio/rom_concepts.hpp>``

.. literalinclude:: ../../../../include/pressio/rom/concepts/fom_steady_with_jac_action.hpp
   :language: cpp
   :lines: 54-115

Semantic requirements
---------------------

Given an instance ``A`` of type ``T`` and an instance ``operand``
of type ``JacobianActionOperandType``,
``SteadyFomWithJacobianAction<T, JacobianActionOperandType>``
is modeled if it is satisfied and all of the following are true:

- methods are blocking: all temporary allocations and operations
  needed to execute those methods
  are completed and no outstanding work remains upon return

- methods only modify non-constant arguments, while const arguments are not modified

- ``auto r = A.createResidual()`` and ``auto result = A.createResultOfJacobianActionOn(operand)``
  return objects with all their "elements" zero initialized

- non-aliasing instantation: this means that doing

  .. code-block:: cpp

     auto r1 = A.createResidual();
     auto r2 = A.createResidual();
     //...
     auto rN = A.createResidual();

  implies that ``r1, r2, ..., rN`` must be distinct objects,
  and such that any modification to ``r1`` does not affect any of the others
  and viceversa.
  In other words, calling ``A.createResidual()`` yields independent instances.
  And the same holds for ``A.createResultOfJacobianActionOn(operand)``.

- ``A.residualAndJacobianAction(state, r, operand, ja, computeJacAction)``

  - overwrites ``r`` with the residual computation

  - overwrites ``ja`` with the result of left-applying the jacobian to ``operand``

  - recomputes the jacobian action only if ``computeJacobian == true``

  - is equality preserving, i.e. equal inputs imply equal outputs

  - let ``J`` represent the Jacobian which we compute the action of,
    then ``J`` must be the jacobian of the residual evaluated for the
    same ``state``. In other words, the Jacobian used for
    computing its action must be mathematically "consistent" with the residual.



Syntax-only example
-------------------

.. code-block:: cpp

   class SampleClass
   {
     public:
       using state_type    = Tpetra::Vector<>; // uses default template parameters
       using residual_type = state_type;

       residual_type createResidual() const;

       Tpetra::MultiVector<> createResultOfJacobianActionOn(const Tpetra::MultiVector<> & operand) const;

       void residual(const state_type & /*state*/,
                     residual_type &    /*result*/) const;

       void residualAndJacobianAction(const state_type & /*state*/,
                                      residual_type & /*residual*/
		                      const Tpetra::MultiVector<> & /*operand*/,
                                      const Tpetra::MultiVector<> & /*result*/,
				      bool /*computeJacAction*/) const;
   }


.. Assuming the default scalar type of Tpetra is ``double``,
   the class above satisfies: ``static_assert(pressio::rom::SteadyFomWithJacobianAction<SampleClass,
   Tpetra::MultiVector<>>, "");``
