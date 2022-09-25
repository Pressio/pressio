.. include:: ../../mydefs.rst

``SemiDiscreteFomWithJacobianAction``
=====================================

Header: ``<pressio/rom_concepts.hpp>``

Namespace: ``pressio::rom``

.. literalinclude:: ../../../../include/pressio/rom/concepts/semi_discrete_fom_with_jac_action.hpp
   :language: cpp
   :lines: 58-80

Semantic requirements
---------------------

:red:`finish`

..
   The concept is modeled only if it is satisfied,
   all subsumed concepts are modeled and given an instance ``S``
   of the trial subspace, and an instance ``A`` of type ``T``,
   then all of the following hold:

   - *non aliasing instantiation*: given the following:

     .. code-block:: cpp

	auto & basis = S.basisOfTranslatedSpace();
	auto ja1 = A.createApplyJacobianResult(basis);
	auto ja2 = A.createApplyJacobianResult(basis);

     ``ja1`` and ``ja2`` must be distinct objects, ``std::addressof(ja1) != std::addressof(ja2)``,
     and such that any modification to ``ja1`` does not affect ``ja2``

   - *blocking operations*: all methods are blocking

   - *equality preserving*

   - *const correctness*: methods may modify only the non-constant operands.
     Operands that are constant must not be modified.

   - a residual instance and the result of the Jaobian action
     must be dimensionally consistent
