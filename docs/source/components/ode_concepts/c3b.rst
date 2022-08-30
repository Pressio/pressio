
``ConstantMassMatrixOperator``
==============================

.. code-block:: cpp

   struct SyntaxOnly
   {
     using mass_matrix_type = /* your type */;

     mass_matrix_type createMassMatrix() const;
     void massMatrix(mass_matrix_type &) const;
   };
