
// this is here so that we have a single place where
// we put the syntax and use it easily in other places
// without duplicating things

struct SyntaxOnly
{
  using time_type            = /* your type */;
  using state_type           = /* your type */;
  using right_hand_side_type = /* your type */;

  state_type createState() const;
  right_hand_side_type createRightHandSide() const;

  void rightHandSide(const state_type &,
		     const time_type &,
		     right_hand_side_type &) const;
};

struct SyntaxOnly
{
  using time_type            = /* your type */;
  using state_type           = /* your type */;
  using right_hand_side_type = /* your type */;

  state_type createState() const;
  right_hand_side_type createRightHandSide() const;

  void rightHandSide(const state_type &,
                     const time_type &,
                     right_hand_side_type &) const;

  template<class BasisType>
  /*return type*/ createApplyJacobianResult(const BasisType & ) const;

  template<class BasisType, class ResultType>
  void applyJacobian(const state_type &,
		     const BasisType &,
		     const time_type &,
		     ResultType &) const;
};

class SyntaxOnly
{
  public:
    using state_type    = /* your type */;
    using residual_type = /* your type */;

  public:
    residual_type createResidual() const;

    template<class BasisType>
    /*return type*/ createApplyJacobianResult(const BasisType & B) const;

    void residual(const state_type &,
		  residual_type &) const;

    template<class BasisType, class ResultType>
    void applyJacobian(const state_type & state,
		       const BasisType &,
		       ResultType &) const;
};


// struct SyntaxOnly
// {
//   using time_type = /* ... */;
//   using state_type                = /* ... */;
//   using discrete_residual_type    = /* ... */;

//   state_type createState() const;
//   discrete_residual_type createDiscreteResidual() const;

//   template<class BasisType>
//   operand_t createResultOfDiscreteTimeJacobianActionOn(const BasisType &) const;

//   // accepting 2 states
//   template<class StepCountType>
//   void discreteResidual(StepCountType,
//                         const time_type & /**/,
//                         const time_type & /**/,
//                         discrete_residual_type & /**/,
//                         const state_type & /**/,
// 			const state_type & /**/) const;

//   template<class StepCountType, class BasisType>
//   void applyDiscreteJacobian(...)

//   // accepting 3 states

//   // accepting 4 states
// };


// struct SyntaxOnly
// {
//   using time_type            = /* your type */;
//   using state_type           = /* your type */;
//   using right_hand_side_type = /* your type */;

//   template<class OperandType>
//   /*return type*/ createApplyMassMatrixResult(const OperandType & ) const;

//   template<class OperandType, class ResultType>
//   void applyMassMatrix(const state_type &,
// 		       const OperandType &,
// 		       const time_type &,
// 		       ResultType & ) const;
// };

// struct SyntaxOnly
// {
//   using time_type            = /* your type */;
//   using state_type           = /* your type */;
//   using right_hand_side_type = /* your type */;

//   state_type createState() const;
//   right_hand_side_type createRightHandSide() const;

//   void rightHandSide(const state_type &,
//                      const time_type &,
//                      right_hand_side_type &) const;

//   template<class BasisType>
//   /*return type*/ createApplyJacobianResult(const BasisType & ) const;

//   template<class BasisType, class ResultType>
//   void applyJacobian(const state_type &,
// 		     const BasisType &,
// 		     const time_type &,
// 		     ResultType & ) const;

//   template<class OperandType>
//   /*return type*/ createApplyMassMatrixResult(const OperandType & ) const;

//   template<class OperandType, class ResultType>
//   void applyMassMatrix(const state_type &,
// 		       const OperandType &,
// 		       const time_type &,
// 		       ResultType & ) const;
// };
