
// this is here so that we have a single place where
// we put the syntax and use it easily in other places
// without duplicating things

struct SyntaxOnly
{
  using residual_operand_type = /* ... */;
  using jacobian_action_operand_type  = /* ... */;

  template<class ResultType>
  void operator()(const residual_operand_type &, ResultType &) const;

  template<class ResultType>
  void operator()(const jacobian_action_operand_type &, ResultType &) const;
}

struct SyntaxOnly
{
  using time_type = /* ... */;
  using right_hand_side_operand_type = /* ... */;

  template<class ResultType>
  void operator()(const right_hand_side_operand_type &, const time_type &, ResultType &) const;
}

struct SyntaxOnly
{
  using time_type = /* ... */;
  using right_hand_side_operand_type = /* ... */;
  using jacobian_action_operand_type  = /* ... */;

  template<class ResultType>
  void operator()(const right_hand_side_operand_type &, const time_type &, ResultType &) const;

  template<class ResultType>
  void operator()(const jacobian_action_operand_type &, const time_type &, ResultType &) const;
}
