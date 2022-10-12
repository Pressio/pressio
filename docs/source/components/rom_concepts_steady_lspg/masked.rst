.. include:: ../../mydefs.rst

``ComposableIntoMaskedProblem``
===============================

Header: ``<pressio/rom_concepts.hpp>``

.. important::

   This concept belongs to the namespace: ``pressio::rom::lspg::steady``

.. literalinclude:: ../../../../include/pressio/rom/concepts/lspg_steady_masked.hpp
   :language: cpp
   :lines: 54-78


Semantic requirements
---------------------

:red:`finish`


.. - if ``N`` is the number of degrees of freeem of the FOM,
..   applying the masker yields operators of size ``n < N``, i.e. the masker
..   in practice subselects elements of the FOM operators
