.. role:: raw-html-m2r(raw)
   :format: html

Steady default problem
======================

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
     class FomReferenceStateType
     >
   auto create_default_steady_problem(const FomSystemType & fomSystem,
                                      DecoderType & decoder,
                                      const RomStateType & romState,
                                      const FomReferenceStateType & fomRefState);

   template<
     class FomSystemType,
     class DecoderType,
     class RomStateType,
     class FomReferenceStateType,
     class PreconditionerType
     >
   auto create_default_steady_problem(const FomSystemType & fomSystem,
                                      DecoderType & decoder,
                                      const RomStateType & romState,
                                      const FomReferenceStateType & fomRefState,
                                      const PreconditionerType & preconditioner);

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

Solve
-----

.. code-block:: cpp

   int main()
   {
   // we assume the rom_state, decoder, fom_system, fom_reference_state already exist

   namespace plspg = pressio::rom::lspg;

   auto problem = plspg::create_default_steady_problem(fom_system, decoder,
                                                       rom_state, fom_reference_state);

   // create nonlinear least-squares solver, for example:
   auto nonlinsolver = pressio::ode::create_gauss_newton(problem, rom_state, ...);
   nonlinsolver.solve(system, rom_state);
   //...
   }
