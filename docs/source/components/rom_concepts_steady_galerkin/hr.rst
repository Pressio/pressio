.. include:: ../../mydefs.rst

``ComposableIntoHyperReducedProblem``
=====================================

Header: ``<pressio/rom_concepts.hpp>``

.. important::

   To avoid name conflicts, this belongs to the namespace: ``pressio::rom::galerkin::steady``


.. literalinclude:: ../../../../include/pressio/rom/concepts/galerkin_steady_hyperreduced.hpp
   :language: cpp
   :lines: 58-92

Semantic and complexity requirements
------------------------------------

:red:`finish`

..
   The concept is modeled only if it is satisfied,
   all subsumed concepts are modeled and letting
   ``S`` be an instance of ``TrialSubspaceType``,
   then all of the following hold:

   - a reduced residual instance ``rGal`` has extent equal to the dimensionality of the subspace, i.e.
     ``pressio::ops::extent(rGal, 0) == S.dimension()``

   - a reduced dense jacobian instance ``JGal`` is such that
     ``pressio::ops::extent(JGal, 0) == pressio::ops::extent(JGal, 1) == S.dimension()``
     must be true, i.e. ``JGal`` must be a square matrix of size
     equal to the dimensionality of the subspace

   - the hyper-reducer operates on the same elements of ``rFom`` and ``JaFom``,
     i.e. it operates on the same degrees of freedom

   - the hyper-reducer represents an approximation of the projection
     of the full FOM residual and jacobian action on the basis

   - *blocking operations*: all methods are blocking

   - *equality preserving*: equal inputs yield equal outputs,
     therefore the hyper-reducer must be invariant

   - *const correctness*: methods may modify only the non-constant operands.
     Operands that are constant must not be modified.

   - :red:`finish complexity saying that a hyper-reducer kind of makes sense if it
     operators on a smaller number of degrees of freedom`



..
   Syntax only
   -----------

   Let ``TrialSubspaceType`` the type of a subspace class
   modeling the `PossiblyAffineTrialColumnSubspace concept <c10.html>`__
   and ``using basis_matrix_type = typename TrialSubspaceType::basis_matrix_type``,
   then we have:

   .. literalinclude:: ./syntax_only_fom_system_concepts.cc
      :language: cpp
      :lines: 43-64
