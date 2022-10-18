.. include:: ../../mydefs.rst

``SystemWithRhsJacobianMassMatrix``
===================================

.. literalinclude:: ../../../../include/pressio/ode/concepts/system_rhs_jacobian_massmatrix.hpp
   :language: cpp
   :lines: 52-105


Semantic requirements
---------------------

:red:`finish`


Syntax-only example
-------------------

.. code-block:: cpp

   class SampleClass
   {
     public:
       using independent_variable_type = double;
       using state_type                = Eigen::VectorXd;
       using right_hand_side_type      = Eigen::VectorXd;
       using jacobian_type             = Eigen::SparseMatrix<double>;
       using mass_matrix_type          = Eigen::SparseMatrix<double>;

       state_type createState() const;
       right_hand_side_type createRightHandSide() const;
       jacobian_type createJacobian() const;
       mass_matrix_type createMassMatrix() const;

       void operator()(const state_type &        /*state*/,
                       independent_variable_type /*ivEval*/,
		       right_hand_side_type &    /*rhsToOverwrite*/,
		       mass_matrix_type &        /*mmToOverwrite*/,
		       jacobian_type &           /*jToOverwrite*/,
		       bool			 /*computeJacobian*/) const;
   }
