
// this is here so that we have a single place where
// we put the syntax and use it easily in other places
// without duplicating things

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

struct SyntaxOnly
{
  using mass_matrix_type = /* your type */;

  mass_matrix_type createMassMatrix() const;
  void massMatrix(mass_matrix_type &) const;
};
