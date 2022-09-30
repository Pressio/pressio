.. include:: ../../mydefs.rst

``FullyDiscreteFomWithJacobianAction``
======================================

Header: ``<pressio/rom_concepts.hpp>``

.. literalinclude:: ../../../../include/pressio/rom/concepts/fom_fully_discrete_with_jac_action.hpp
   :language: cpp
   :lines: 58-115

Semantic requirements
---------------------

:red:`finish`


Syntax-only example
-------------------

.. code-block:: cpp

   class SampleClass
   {
     public:
       using time_type              = double;
       using state_type             = Tpetra::Vector<>; // uses default template parameters
       using discrete_residual_type = state_type;

       discrete_residual_type createDiscreteTimeResidual() const;

       Tpetra::MultiVector<>
       createResultOfDiscreteTimeJacobianActionOn(const Tpetra::MultiVector<> & /*operand*/);

       void discreteTimeResidualAndJacobianAction(pressio::ode::StepCount::value_type /*step number*/,
						  time_type /*evalTime*/,
						  time_type /*dt*/,
						  discrete_residual_type & /*residual*/,
						  const Tpetra::MultiVector<> & /*operand*/,
						  bool /*computeJacobianAction*/,
						  Tpetra::MultiVector<> & /*jacobianAction*/,
						  const state_type & /*currentState*/,
						  const state_type & /*previousState*/) const;
   }


Assuming the default scalar type of Tpetra is ``double``,
the class above satisfies: ``static_assert(pressio::rom::FullyDiscreteFomWithJacobianAction<SampleClass, 2,
Tpetra::MultiVector<>>, "");``
