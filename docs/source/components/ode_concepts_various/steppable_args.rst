
``SteppableWithAuxiliaryArgs``
==============================

The concept ``SteppableWithAuxiliaryArgs`` specifies that
a type ``T`` represents something that can perform
a single step and to do so it needs at least one auxiliary object.

Syntax only
-----------

.. literalinclude:: ./syntax_various.cc
   :language: cpp
   :lines: 20-33

..
   Full concept
   ------------

   .. code-block:: cpp

       template <class T, class AuxT, class ...Args>
       concept SteppableWithAuxiliaryArgs =
	 //
	 // purely syntactic requirements
	 //
	 requires(){
	   typename T::state_type;
	   typename T::independent_variable_type;
	   requires pressio::Time<typename T::independent_variable_type>;
	 } &&
	 requires(const T & A,
		  const typename T::state_type & s,
		  pressio::ode::StepStartAt<typename T::independent_variable_type> sst,
		  pressio::ode::StepCount stepNumber,
		  pressio::ode::StepSize<typename T::independent_variable_type> dt,
		  AuxT && arg1,
		  Args && ... args)
	 {
	   A(s, sst, stepNumber, dt,
	     std::forward<AuxT>(arg1),
	     std::forward<Args>(args)...);
	 } &&

	 //
	 // execution/language axioms
	 //
	 axiom TimeAndAuxiliaryArgsDependence(){
	   // state has a meaningful dependence on time
	   // and the auxiliary args are needed to perform a step
	 } &&
	 axiom BlockingOperation(){
	   // callable method is blocking (completes before returning)
	 } &&
	 axiom ConsistentUnits(){
	   // units of time and state are consistent with
	   // their dependency relationship
	 } &&
	 axiom ConstCorrectness(){
	   // const qualification is preserved, methods do NOT modify const arguments
	 };

``StronglySteppableWithAuxiliaryArgs``
======================================

.. code-block:: cpp

    template <class T, class AuxT, class ...Args>
    concept StronglySteppableWithAuxiliaryArgs =
       SteppableWithAuxiliaryArgs<T, AuxT, Args...> &&
       axiom StronglyGuaranteeingStep(){
         // doing one step via operator() is strongly guaranteeing:
	 // if an exception is thrown inside operator(), upon return,
	 // the state object is guaranteed to be the same as
	 // before starting the step
       };
