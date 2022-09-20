
// this is here so that we have a single place where
// we put the syntax and use it easily in other places
// without duplicating things



struct SyntaxOnly
{
  using basis_matrix_type  = /**/;

  const basis_matrix_type & basis() const;
  std::size_t dimension() const;
  bool isColumnSpace() const;
  bool isRowSpace() const;
}

struct SyntaxOnly
{
  using reduced_state_type = ;
  using full_state_type    = ;
  using basis_matrix_type  = ;

  reduced_state_type createReducedState() const;
  full_state_type    createFullState() const;
  void               mapFromReducedState(const reduced_state_type &, full_state_type &) const;
  full_state_type    createFullStateFromReducedState(const reduced_state_type &) const;

  const basis_matrix_type & basisOfTranslatedSpace() const;
  const full_state_type   & translationVector() const;
  const basis_matrix_type & basis() const;
  std::size_t dimension() const;
  bool isColumnSpace() const;
  bool isRowSpace() const;
}



// struct SyntaxOnly
// {
//   using reduced_state_type = ;
//   using full_state_type    = ;
//   using tangent_space_basis_matrix_type = ;
//   using offset_type        = ;

//   reduced_state_type createReducedState() const;
//   full_state_type    createFullState() const;
//   void               mapFromReducedState(const reduced_state_type &, full_state_type &) const;
//   full_state_type    createFullStateFromReducedState(const reduced_state_type &) const;
//   const tangent_space_basis_matrix_type & viewTangentSpaceBasis() const;
//   const offset_type & viewAffineOffset() const;
//   std::size_t dimension() const;
// }
