
concepts: system
================

.. toctree::
   :maxdepth: 1

   ode_concepts_system/rhs
   ode_concepts_system/rhs_jac
   ode_concepts_system/rhs_massmatrix
   ode_concepts_system/rhs_jac_massmatrix
   ode_concepts_system/fully_discrete

concepts: various
=================

.. toctree::
    :maxdepth: 1

    ode_concepts_various/steppable
    ode_concepts_various/steppable_args
    ode_concepts_various/step_size_pol
    ode_concepts_various/step_size_pol_with_scaling
    ode_concepts_various/state_observer
    ode_concepts_various/state_guesser



.. A fundamental part of pressio is that concept-driven design has been
   a key part of the development since its early stages.
   Concepts have formally entered the C++ standard with C++20, but are
   such a new feature that it will require some time before they are more widely adopted.

   We clearly state upfront that we embrace the following definitions:

   - Constraints define the statically evaluable predicates on the properties and syntax
   of types, but do not represent cohesive abstractions.

   - Axioms state semantic requirements on types that should not be statically evaluated.
   An axiom is an invariant that is assumed to hold (as opposed required to be checked)
   for types that meet a concept. Axioms specify the meaning of those interfaces and relationships, type invariants, and complexity guarantees

   – Concepts are predicates that represent general, abstract, and stable requirements of
   generic algorithms on their argument. They are defined in terms of constraints and
   axioms

   **The distinctive property that separates a
   concept from a constraint is that it has semantic properties**
