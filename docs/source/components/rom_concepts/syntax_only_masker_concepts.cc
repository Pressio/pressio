
// this is here so that we have a single place where
// we put the syntax and use it easily in other places
// without duplicating things

struct TimeInvariantMasker
{
  using operand_type = /* ... */;
  using result_type  = /* ... */;

  result_type createApplyMaskResult(const operand_type &) const;
  void operator()(const operand_type &, result_type &) const;
}

struct TimeDependentMasker
{
  using operand_type = /* ... */;
  using result_type  = /* ... */;

  result_type createApplyMaskResult(const operand_type &) const;

  template<class TimeType>
  void operator()(const operand_type &, const TimeType &, result_type &) const;
}
