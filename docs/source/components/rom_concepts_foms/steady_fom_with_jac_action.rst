.. include:: ../../mydefs.rst

``SteadyFomWithJacobianAction``
===============================

Header: ``<pressio/rom_concepts.hpp>``

Namespace: ``pressio::rom``

.. literalinclude:: ../../../../include/pressio/rom/concepts/steady_fom_with_jac_action.hpp
   :language: cpp
   :lines: 56-77

Semantic requirements
---------------------

:red:`finish`

..
   The concept is modeled only if it is satisfied,
   all subsumed concepts are modeled and all of the following hold:

   - *non aliasing instantiation*: given the following:

     .. code-block:: cpp

	auto r1 = A.createResidual();
	auto r2 = A.createResidual();

     ``r1`` and ``r2`` must be distinct objects, ``std::addressof(r1) != std::addressof(r2)``,
     and such that any modification to ``r1`` does not affect ``r2``

   - *blocking operations*: all methods are blocking

   - *equality preserving*

   - *const correctness*: methods may modify only the non-constant operands.
     Operands that are constant must not be modified.

   - a residual instance and the result of the Jaobian action
     must be dimensionally consistent
