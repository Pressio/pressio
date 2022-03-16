.. role:: raw-html-m2r(raw)
   :format: html

Unsteady default problem
========================

.. note::

    Defined in: ``<pressio/rom_lspg.hpp>``

    Public namespace: ``pressio::rom::lspg``

.. warning::

    Prerequisite reading:
    Before you read this page, make sure you
    read the `overall design idea of the unsteady LSPG <rom_lspg_unsteady.html>`_.

API
---

.. code-block:: cpp

   // overload for continuous-time systems
   template<
     class FomSystemType,
     class DecoderType,
     class RomStateType,                                                             (1)
     class FomReferenceStateType
     >
   ReturnType create_default_unsteady_problem(pressio::ode::StepScheme scheme,
                                              const FomSystemType & fomSystem,
                                              DecoderType & decoder,
                                              const RomStateType & romState,
                                              const FomReferenceStateType & fomRefState);

   // overload for discrete-time systems
   template<
     std::size_t num_states,
     class FomSystemType,
     class DecoderType,
     class RomStateType,                                                             (2)
     class FomReferenceStateType
     >
   ReturnType create_default_unsteady_problem(const FomSystemType & fomSystem,
                                              DecoderType & decoder,
                                              const RomStateType & romState,
                                              const FomReferenceStateType & fomRefState);

Parameters and Requirements
^^^^^^^^^^^^^^^^^^^^^^^^^^^

* 
  ``fomSystem``\ :

  * instance of your FOM adapter specifying the FOM problem
  * for 1: must satisfy the `continuous-time API <rom_fom_apis.html>`_
  * for 2: must satisfy the `discrete-time API <rom_fom_apis.html>`_

* 
  ``decoder``\ :

  * decoder object
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

* ``num_states``\ :

  * only needed for the discrete-time case
  * *total* number of states you need to use (must be <= 3), if you need more,
    please tell us by `openining an issue <https://github.com/Pressio/pressio/issues>`_
