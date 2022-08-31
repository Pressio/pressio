
``StepSizePolicyWithScalingFactor``
===================================

This concept specifies that a type ``T`` represents a policy
to set the step size
to use at a given step as well as the logic on how to
potentially "reduce/scale" the step size for a given step.

Syntax only
-----------

.. literalinclude:: ./syntax_various.cc
   :language: cpp
   :lines: 46-55

..
   Full concept
   ------------

   .. code-block:: cpp

       template <class T>
       concept StepSizePolicyWithReductionScheme =
	 //
	 // purely syntactic requirements
	 //
	 requires(const T & A,
		  pressio::ode::StepCount stepNumber,
		  const pressio::ode::StepStartAt<typename T::time_type> & sst,
		  pressio::ode::StepSize<IndepVarType> & dt,
		  pressio::ode::StepSizeMinAllowedValue<IndepVarType> & dtMin,
		  pressio::ode::StepSizeScalingFactor<IndepVarType> & factor)
	 {
	   A(stepNumber, sst, dt, dtMin, factor);
	 } &&

	 //
	 // execution/language axioms
	 //
	 axiom BlockingOperation(){
	   // callable method is blocking (completes before returning)
	 } &&
	 axiom ConsistentUnits(){
	   // step units are consistent with time
	 } &&
	 axiom SOMETHINGELSE(){

	 } &&
	 axiom ConstCorrectness(){
	   // const qualification is preserved, methods do NOT modify const arguments
	 };
