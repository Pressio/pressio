.. role:: raw-html-m2r(raw)
   :format: html

Iteratively reweighted least squares
====================================

Defined in header ``<pressio/solvers_nonlinear.hpp>``


API, Parameters and Requirements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: cpp

   namespace pressio{ namespace nonlinearsolvers{ namespace experimental{

   template<class ProblemClassType, class LinearSolverType>
   auto create_irls_gauss_newton(const ProblemClassType & system,
		                 LinearSolverType && lsolver);

   }}}

* ``system``:

  - instance of your problem class defining the problem

  - must satisfy either the `SystemWithResidualAndJacobian <nonlinearsolvers_concepts/c1.html>`__,
    or `SystemWithFusedResidualAndJacobian <nonlinearsolvers_concepts/c2.html>`__

* ``lsolver``:

  * linear solver for solving the normal equations, choose one from `linear solvers <linsolvers.html>`_
  * if you want to implement your own, then it has to conform to the `linear solver API <linsolvers.html>`_


.. warning::

   Note that this is still inside the experimental namespace.


Example usage
^^^^^^^^^^^^^

something
