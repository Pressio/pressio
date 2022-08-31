
``StateGuesser``
================

This concept specifies that a type ``T``
represents a valid guesser called to guess the state *before*
a step is performed.

Syntax only
-----------

.. literalinclude:: ./syntax_various.cc
   :language: cpp
   :lines: 68-75

..
   Full concept
   ------------

   .. code-block:: cpp

       template <class T>
       concept StateGuesser =
	 //
	 // purely syntactic requirements
	 //
	 template<class IndepVarType, class StateType>
	 requires(const T & A,
		  pressio::ode::StepCount stepNumber,
		  pressio::ode::TimeValue<IndepVarType> currTime,
		  StateType & s)
	 {
	   A(stepNumber, currTime, s);
	 } &&

	 //
	 // execution/language axioms
	 //
	 axiom BlockingOperation(){
	   // callable method is blocking (completes before returning)
	 } &&
	 axiom ConstCorrectness(){
	   // const qualification is preserved, methods do NOT modify const arguments
	 } &&
	 axiom ConsistentModification(){
	   // state is overwritten such that the new state
	   // is consistent with the previous one
	 };
