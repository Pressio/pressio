
``HyperReduceableAndMaskableWith``
==================================

Header: ``<pressio/rom_concepts.hpp>``

Namespace: ``pressio::rom::galerkin::steady``

.. literalinclude:: ../../../../include/pressio/rom/concepts/steady_galerkin_hyperreduceable_maskable_with.hpp
   :language: cpp
   :lines: 59-108

Semantic and complexity requirements
------------------------------------

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

- the hyper-reducer touches the same "elements" of the masked residual and
  jacobian action, i.e. it operates on the same degrees of freedom

- the hyper-reducer implements an approximation of the projection
  of the *masked* FOM residual and jacobian action on the basis

- *blocking operations*: all methods are blocking

- *equality preserving*: therefore it is invariant

- *const correctness*: methods may modify only the non-constant operands.
  Operands that are constant must not be modified.

- if ``N`` is the number of degrees of freeem of the FOM,
  applying the masker yields operators of size ``n < N``, i.e. the masker
  in practice subselects elements of the FOM operators

- :red:`finish discussing the complexity argument`


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
