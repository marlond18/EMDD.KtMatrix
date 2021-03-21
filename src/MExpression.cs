using System;
using EMDD.KtNumerics;
using EMDD.KtExpressions;
using System.Collections.Generic;

namespace EMDD.KtMatrix
{
    public class MExpression : MElement
    {
        internal readonly Expression _value;

        public MExpression(Expression expression)
        {
            if (expression == null) throw new ArgumentNullException(nameof(expression));
            _value = expression.Clone();
        }

        public override bool IsZero => _value == 0;

        public override MElement AddTo(MElement other) => other switch
        {
            MNumeric numeric => new MExpression(_value + numeric._value),
            MExpression expression => new MExpression(_value + expression._value),
            _ => other + this,
        };

        public override MElement Negate() => new MExpression(-_value);

        public override MElement MultiplyTo(MElement other) => other switch
        {
            MNumeric numeric => new MExpression(_value * numeric._value),
            MExpression expression => new MExpression(_value * expression._value),
            _ => other * this,
        };

        public IEnumerable<(double location,Number result)> EvaluateOnEqualInterval(int division) => _value.EvaluateOnEqualInterval(division);

        public override MElement Inverse() => new MExpression(1 / _value);
        public static implicit operator Expression(MExpression expression) => expression._value;
        public override string ToString() => _value.ToString();
        public string ToStringPieceWise() => _value.ToStringPieceWise();
        public override int GetHashCode() => _value.GetHashCode() * 16;

        public override bool Equals(MElement other) => other switch
        {
            MExpression expression => _value == expression._value,
            MNumeric num => _value == num._value,
            _ => false
        };

        public override MElement Clone() => new MExpression(_value.Clone());

        public override double ToDouble() => 0;
    }
}