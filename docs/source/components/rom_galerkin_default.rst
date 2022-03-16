.. role:: raw-html-m2r(raw)
   :format: html

Default problem
===============

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
     class FomReferenceStateType
     >                                                                              (1)
   auto create_default_explicit_problem(pressio::ode::StepScheme scheme,
                                        const FomSystemType & fomSystem,
                                        DecoderType & decoder,
                                        const RomStateType & romState,
                                        const FomReferenceStateType & fomRefState)

   template<
     class FomSystemType,
     class DecoderType,
     class RomStateType,
     class FomReferenceStateType
     >                                                                              (2)
   auto create_default_implicit_problem(pressio::ode::StepScheme scheme,
                                        const FomSystemType & fomSystem,
                                        DecoderType & decoder,
                                        const RomStateType & romState,
                                        const FomReferenceStateType & fomRefState)

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

  * instance of your FOM adapter specifying the FOM problem
  * Must satisfy one of the APIs suitable for Galerkin, see `API list <rom_fom_apis.html>`_

* 
  ``decoder``\ :

  * your decoder object
  * must satify the requirements listed `here <rom_decoder.html>`_

* 
  ``romState``\ :

  * currently, it must be either an Eigen vector or a Kokkos 1D view

* 
  ``fomRefState``\ :

  * your FOM reference state that is used when reconstructing the FOM state
  * must be copy-constructible and the following must be true:

    .. code-block:: cpp

       std::is_same<FomReferenceStateType, typename DecoderType::fom_state_type>::value == true

Galerkin sample code:
^^^^^^^^^^^^^^^^^^^^^

.. tip::

    - `Explicit Galerkin sample code 1 <https://pressio.github.io/pressio-tutorials/known_types/proms_galerkin_default_explicit.html>`_
    - `Implicit Galerkin sample code 1 <https://pressio.github.io/pressio-tutorials/known_types/proms_galerkin_default_implicit.html>`_
