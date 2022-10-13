.. role:: raw-html-m2r(raw)
   :format: html

Gauss-Newton (via QR factorization)
===================================

Defined in header ``<pressio/solvers_nonlinear.hpp>``


API, Parameters and Requirements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: cpp

   namespace pressio{ namespace nonlinearsolvers{

   template<class ProblemClassType, class QRSolverType>
   auto create_gauss_newton(const ProblemClassType & system,
                            QRSolverType && qrsolver);

   }}

* ``system``

  - instance of your problem class defining the problem

  * accepts a system satisfying either the
    `SystemWithResidualAndJacobian <nonlinearsolvers_concepts/c1.html>`__ or
    `SystemWithFusedResidualAndJacobian <nonlinearsolvers_concepts/c2.html>`__

* ``qrsolver``

  * solver needed to solve the QR-based formulation of the least-squares problem `see this <https://en.wikipedia.org/wiki/QR_decomposition>`_
  * we suggest to use the `pressio QR package <qr.html>`_
  * if you want to implement your own, then it has to conform to the `this API <qr.html>`_


Example usage
^^^^^^^^^^^^^

bla
