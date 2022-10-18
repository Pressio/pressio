.. include:: ../../mydefs.rst

``SystemWithRhsAndJacobian``
============================

.. literalinclude:: ../../../../include/pressio/ode/concepts/system_rhs_jacobian.hpp
   :language: cpp
   :lines: 52-100


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

       state_type createState() const;
       right_hand_side_type createRightHandSide() const;
       jacobian_type createJacobian() const;

       void operator()(const state_type &        /*state*/,
                       independent_variable_type /*ivEval*/,
		       right_hand_side_type &    /*rhsToOverwrite*/,
		       jacobian_type &           /*jToOverwrite*/,
		       bool			 /*computeJacobian*/) const;
   }
