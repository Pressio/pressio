
``MassMatrixOperator``
======================

.. code-block:: cpp

   struct SyntaxOnly
   {
     using independent_variable_type = /* your type */;
     using state_type                = /* your type */;
     using mass_matrix_type          = /* your type */;

     mass_matrix_type createMassMatrix() const;

     void massMatrix(const state_type &,
		     const independent_variable_type &,
		     mass_matrix_type &) const;
   };
