using EMDD.KtNumerics;

namespace EMDD.KtMatrix
{
    public class MNumeric : MElement
    {
        internal readonly Number _value;

        public MNumeric(Number value)
        {
            _value = value;
        }

        public override bool IsZero => _value == 0;
        public override MElement AddTo(MElement other) => other is MNumeric number ? _value + number._value : other + this;
        public override MElement Clone() => new MNumeric(_value);
        public override bool Equals(MElement other) => other is MNumeric numeric && _value == numeric._value;
        public override int GetHashCode() => _value.GetHashCode() * 32;
        public override MElement Inverse() => new MNumeric(1 / _value);
        public override MElement MultiplyTo(MElement other) => other is MNumeric numeric ? new MNumeric(_value * numeric._value) : other * this;
        public override MElement Negate() => new MNumeric(-_value);
        public override double ToDouble() => _value.ToDouble();
        public override string ToString() => _value.ToString();
        public static implicit operator Number(MNumeric num) => num._value;
    }
}
