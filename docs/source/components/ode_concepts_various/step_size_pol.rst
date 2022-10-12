.. include:: ../../mydefs.rst

``StepSizePolicy``
==================

.. literalinclude:: ../../../../include/pressio/ode/concepts/ode_step_size_policy.hpp
   :language: cpp
   :lines: 52-66

Semantic requirements
---------------------

:red:`finish`

..
   axiom BlockingOperation(){
     // callable method is blocking (completes before returning)
   } &&
   axiom ConsistentUnits(){
     // step units are consistent with time
   } &&
   axiom ConstCorrectness(){
     // const qualification is preserved, methods do NOT modify const arguments
   };
