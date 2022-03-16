.. role:: raw-html-m2r(raw)
   :format: html

Hyper-reduced problem
=====================

.. note::

    Defined in: ``<pressio/rom_galerkin.hpp>``

    Public namespace: ``pressio::rom::galerkin``

.. warning::

    Prerequisite reading:
    Before you read this page, make sure you
    read the `overall design idea of Galerkin <rom_galerkin.html>`_.

API
---

.. code-block:: cpp

   template<
     class FomSystemType,
     class DecoderType,
     class RomStateType,
     class FomReferenceStateType,
     class ProjectorType
     >                                                                              (1)
   auto create_hyperreduced_explicit_problem(pressio::ode::StepScheme scheme,
                                             const FomSystemType & fomSystem,
                                             DecoderType & decoder,
                                             const RomStateType & romState,
                                             const FomReferenceStateType & fomRefState,
                                             const ProjectorType & projector);

   template<
     class FomSystemType,
     class DecoderType,
     class RomStateType,
     class FomReferenceStateType,
     class ProjectorType
     >                                                                              (2)
   auto create_hyperreduced_implicit_problem(pressio::ode::StepScheme scheme,
                                             const FomSystemType & fomSystem,
                                             DecoderType & decoder,
                                             const RomStateType & romState,
                                             const FomReferenceStateType & fomRefState,
                                             const ProjectorType & projector);

This function returns an instance of the desired Galerkin problem.

Parameters and Requirements
^^^^^^^^^^^^^^^^^^^^^^^^^^^

* 
  ``scheme``\ :


  * enum value to specify the stepper scheme
  * (1) explicit Galerkin, see `valid enum scheme choices <ode_steppers_explicit.html>`_
  * (2) implicit Galerkin, see `valid enum scheme choices <ode_steppers_implicit.html>`_

* 
  ``fomSystem``\ :

  * instance of your FOM adapter specifying the FOM problem :raw-html-m2r:`<br/>`
  * must satisfy one of the APIs suitable for Galerkin, see `API list <rom_fom_apis.html>`_

* 
  ``decoder``\ :

  * decoder object
  * must satify the requirements listed `here <rom_decoder.html>`_

* 
  ``romState``\ :

  * ROM state
  * currently, it must be either an Eigen vector or a Kokkos 1D view

* 
  ``fomRefState``\ :

  * your FOM reference state that is used when reconstructing the FOM state
  * must be copy-constructible and the following must be true:

    .. code-block:: cpp

       std::is_same<FomReferenceStateType, typename DecoderType::fom_state_type>::value == true

* 
  ``projector``\ :

  * operator responsible for projectng the FOM operators onto the reduced space
  * must be a functor with a specific API, see details below

Projector
---------

\todo: explain what projector is for since it is critical.

The projector must be a functor as follows:

.. code-block:: cpp

   struct ValidProjector
   {
     template<class operand_type, class TimeType, class ResultType>
     void operator()(const operand_type & operand,
                     const TimeType time,
                     ResultType & result) const
     {
       // apply the projection to the operand, store into result
     }
   };

Note that here the ``operator()`` is templated on the operand
because the projector can be called with different arguments depending
on whether you are selected an explicit or implicit Galerkin problem:

* for explicit: the projector *always* receives as operand argument
  an instance of your FOM velocity
* for implicit: the operand is an instance of your FOM velocity
  as well as an instance of the decoder's Jacobian.

You are responsible to handle the various cases.
Obviously, you are not required to template the operand.
If you are working with specific types, you can simply specialize it.
For example, suppose that you are working with an application
that uses for the velocity a type named ``FomVelocityType``\ ,
and for the decoder's Jacobian you are using a type named ``ADenseMatrixType``.
In such case, a valid projector can be written as:

.. code-block:: cpp

   struct ValidProjector
   {
     template<class TimeType, class ResultType>
     void operator()(const FomVelocityType & operand,
                     const TimeType time,
                     ResultType & result) const
     {
       // apply the projection to the operand, store into result
     }

     template<class TimeType, class ResultType>
     void operator()(const ADenseMatrixType & operand,
                     const TimeType time,
                     ResultType & result) const
     {
       // apply the projection to the operand, store into result
     }
   };

The compiler will select the best match.
Note that the result is NOT templated.
In general, you can safely assume that ``result`` is indexable as ``(i,...)``.
\todo describe this more
