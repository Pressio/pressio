
``StateObserver``
=================

This concept specifies that a type ``T``
represents a valid observer to monitor the time-evolving state.

Syntax only
-----------

.. literalinclude:: ./syntax_various.cc
   :language: cpp
   :lines: 58-65

..
   Full concept
   ------------

   .. code-block:: cpp

       template <class T>
       concept StateObserver =
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
	 };
