.. include:: ../../mydefs.rst

``SteppableWithAuxiliaryArgs``
==============================

.. literalinclude:: ../../../../include/pressio/ode/concepts/ode_steppable_with_args.hpp
   :language: cpp
   :lines: 56-77


``StronglySteppableWithAuxiliaryArgs``
======================================

.. literalinclude:: ../../../../include/pressio/ode/concepts/ode_steppable_with_args.hpp
   :language: cpp
   :lines: 79-85

Semantic requirements
---------------------

:red:`finish`

..
  // doing one step via operator() is strongly guaranteeing:
  // if an exception is thrown inside operator(), upon return,
  // the state object is guaranteed to be the same as
  // before starting the step
