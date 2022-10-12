.. include:: ../../mydefs.rst

``Steppable``
=============

.. literalinclude:: ../../../../include/pressio/ode/concepts/ode_steppable.hpp
   :language: cpp
   :lines: 52-71

Semantic requirements
---------------------

..
   axiom ConstCorrectness(){
     // const qualification is preserved, methods do NOT modify const arguments
   } &&
   //
   // stepCount is supposed to be 1,2,3,4,..., so StepCount \in Z^+
   // this is important

:red:`finish`


``StronglySteppable``
=====================

.. literalinclude:: ../../../../include/pressio/ode/concepts/ode_steppable.hpp
   :language: cpp
   :lines: 73-80

Semantic requirements
---------------------

:red:`finish`

..
  // doing one step via operator() is strongly guaranteeing:
  // if an exception is thrown inside operator(), upon return,
  // the state object is guaranteed to be the same as
  // before starting the step
