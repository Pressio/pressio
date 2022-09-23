
``MaskableWith``
================

Header: ``<pressio/rom_concepts.hpp>``

Namespace: ``pressio::rom::lspg::steady``

.. literalinclude:: ../../../../include/pressio/rom/concepts/steady_lspg_maskable_with.hpp
   :language: cpp
   :lines: 59-77

Semantic requirements
---------------------

Let ``S`` be an instance of ``TrialSubspaceType``, then:

- *blocking operations*: all methods are blocking

- *equality preserving*: therefore it is invariant

- *const correctness*: methods may modify only the non-constant operands.
  Operands that are constant must not be modified.

- if ``N`` is the number of degrees of freeem of the FOM,
  applying the masker yields operators of size ``n < N``, i.e. the masker
  in practice subselects elements of the FOM operators

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
