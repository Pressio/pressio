.. role:: raw-html-m2r(raw)
   :format: html

Steady masked problem
=====================

.. note::

    Defined in: ``<pressio/rom_lspg.hpp>``

    Public namespace: ``pressio::rom::lspg``

.. warning::

    Prerequisite reading:
    Before you read this page, make sure you
    read the `overall design idea behind steady LSPG <rom_lspg_steady.html>`_.

API
---

.. code-block:: cpp

   template<
     class FomSystemType,
     class DecoderType,
     class RomStateType,
     class FomReferenceStateType,
     class MaskerType
     >
   auto create_masked_steady_problem(const FomSystemType & fomSystem,
                                     DecoderType & decoder,
                                     const RomStateType & romState,
                                     const FomReferenceStateType & fomRefState,
                                     const MaskerType & masker);

   template<
     class FomSystemType,
     class DecoderType,
     class RomStateType,
     class FomReferenceStateType,
     class MaskerType,
     class PreconditionerType
     >
   auto create_masked_steady_problem(const FomSystemType & fomSystem,
                                     DecoderType & decoder,
                                     const RomStateType & romState,
                                     const FomReferenceStateType & fomRefState,
                                     const PreconditionerType & preconditioner,
                                     const MaskerType & masker);

Parameters and Requirements
^^^^^^^^^^^^^^^^^^^^^^^^^^^

* 
  ``fomSystem``\ :

  * instance of your FOM adapter specifying the FOM problem
  * must satisfy the steady API, see `here <rom_fom_apis.html>`_

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
