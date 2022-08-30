
// this is here so that we have a single place where
// we put the syntax and use it easily in other places
// without duplicating things

struct SyntaxOnly
{
  using reduced_state_type = ;
  using basis_type         = ;
  using full_state_type    = ;

  reduced_state_type createReducedState() const;
  full_state_type createFullState() const;
  void mapFromReducedState(const redStateType &, full_state_type &) const;
  full_state_type createFullStateFromReducedState(const redStateType &) const;
  const basis_type & viewBasis() const;
}

struct SyntaxOnly
{
  using reduced_state_type = ;
  using basis_type         = ;
  using full_state_type    = ;

  reduced_state_type createReducedState() const;
  full_state_type createFullState() const;
  void mapFromReducedState(const redStateType &, full_state_type &) const;
  full_state_type createFullStateFromReducedState(const redStateType &) const;
  const basis_type & viewBasis() const;
  const full_state_type & viewAffineOffset() const;
}
