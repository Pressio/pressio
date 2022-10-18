.. role:: raw-html-m2r(raw)
   :format: html

.. include:: ../mydefs.rst

Policy-based Customization
===========================

This page describes the "policy-based use" of an "implicit stepper"
for systems expressed either as

.. math::
  :label: ode_implicit_system_pol1

    \frac{d \boldsymbol{y}}{dt} =
    \boldsymbol{f}(\boldsymbol{y},t; ...),  \qquad y(t_0) = y_0

or with a mass matrix:

.. math::
   :label: ode_implicit_system_pol2

    M(\boldsymbol{y}, t, ...) \frac{d \boldsymbol{y}}{dt} =
    \boldsymbol{f}(\boldsymbol{y},t; ...),  \qquad y(t_0) = y_0


What do we mean by `policy-based` use?
--------------------------------------

Rather than specifying the problem via a system class as described in the `standard use <ode_steppers_implicit_semidiscrete_standard_use.html>`__, you define your problem using a *policy-based approach*, i.e., you provide pressio with one policy to compute the residual and one for computing the Jacobian.


Why is this useful?
-------------------

This is useful if, given an implicit scheme, you want (or need) to have full control on how you compute the discrete residual and jacobian. This is not just about how the operators are computed: for some problems, e.g. LSPG-ROM, one has to customize to definition of residual and jacobian. :red:`say more here`.

What is the downside?
---------------------

The main downside of this functionality is that you are **fully responsible** of computing things correctly according to the the chosen scheme. Therefore, you need to know the details of scheme. So you have to be careful in doing this.


Usage
-----

:red:`TODO`



..
   This use case practically involves these steps:

   1. you define your problem by defining custom *policy* classes
   2. you use the policies to create a pressio stepper
   3. you use the stepper to advance in time


   1. Define your problem
   ~~~~~~~~~~~~~~~~~~~~~~

   Here you define your problem by providing two classes
   that know how to compute the discrete-time residual and
   the discrete-time jacobian of


   * Your residual policy class must conform to the ``ResidualPolicy`` concept:

   .. literalinclude:: ./ode_concepts/syntax_only_for_all_concepts.cc
      :language: cpp
      :lines: 240-258


   * Your jacobian policy class must conform to ``JacobianPolicy`` concept:

   .. literalinclude:: ./ode_concepts/syntax_only_for_all_concepts.cc
      :language: cpp
      :lines: 260-276


   If you want to use custom policies for computing residual and Jacobian,
   **you are responsible** for ensuring things are correct.
   In particular, you should be aware of the following:

   * ``state``

     * passed to ``call`` operator of the policies, contains the prediction at ``n+1``.

   * ``stencilStates``\ , ``auxRhs``

     * the types of these you don't need to know
     * contain the needed auxiliary states and RHS evaluations, respectively, needed to compute the operators for a certain scheme. All you need to know about these containers is the following:

   .. list-table::
      :widths: 20 80
      :header-rows: 1

      * - Scheme
	- Description/Info
      * - ``BDF1``
	- ``stencilStates``: contains state at n-th step :raw-html-m2r:`<br/>` Usage: ``const auto & yn = stencilStates(pressio::ode::n());`` :raw-html-m2r:`<br/>` :raw-html-m2r:`<br/>` ``auxRhs`` : Empty
      * - ``BDF2``
	- ``stencilStates``: contains states at n-th and (n-1)-th step :raw-html-m2r:`<br/>` Usage: ``const auto & yn = stencilStates(pressio::ode::n());`` :raw-html-m2r:`<br/>` :raw-html-m2r:`&nbsp; &nbsp; &nbsp; &nbsp; &emsp` ``const auto & ynm1 = stencilStates(pressio::ode::nMinusOne());`` :raw-html-m2r:`<br/>` :raw-html-m2r:`<br/>` ``auxRhs`` : Empty
      * - ``CrankNicolson``
	- ``stencilStates``: contains states at n-th step :raw-html-m2r:`<br/>` Usage: ``const auto & yn = stencilStates(pressio::ode::n());`` :raw-html-m2r:`<br/>` :raw-html-m2r:`<br/>` ``auxRhs``: contains evaluations of the RHS are n-th and (n+1)-th steps :raw-html-m2r:`<br/>` Usage: ``auto & fn = auxRhs(pressio::ode::n());`` :raw-html-m2r:`<br/>` :raw-html-m2r:`&nbsp; &nbsp; &nbsp; &nbsp; &emsp` ``auto & fnp1 = auxRhs(pressio::ode::nPlusOne());``



   2. Instantiate a stepper
   ~~~~~~~~~~~~~~~~~~~~~~~~

   .. code-block:: cpp

       // overload 2
       template<class ResidualPolicyType, class JacobianPolicyType>
       auto create_implicit_stepper(pressio::ode::StepScheme name,
				    ResidualPolicyType && residualPolicy,
				    JacobianPolicyType && jacobianPolicy);

   *
     ``scheme_name``

     * one of the following enum values ``pressio::ode::StepScheme::{BDF1, BDF2, CrankNicolson}``.


   3. Use the stepper
   ~~~~~~~~~~~~~~~~~~

   Using the stepper to solve this problem can be done similarly
   to example shown above in the section at the top.
