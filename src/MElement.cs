using EMDD.KtNumerics;
using EMDD.KtExpressions;

namespace EMDD.KtMatrix
{
    public abstract class MElement
    {
        public abstract MElement Clone();
        public abstract bool IsZero { get; }
        public abstract MElement AddTo(MElement other);
        public abstract MElement Negate();
        public abstract MElement MultiplyTo(MElement other);
        public abstract MElement Inverse();
        public abstract override string ToString();
        public abstract override int GetHashCode();
        public abstract bool Equals(MElement other);
        public abstract double ToDouble();

        public static implicit operator MElement(double val) => (Number)val;
        public static implicit operator MElement(Number val) => new MNumeric(val);
        public static implicit operator MElement(Expression val) => new MExpression(val);

        public override bool Equals(object obj)
        {
            if (ReferenceEquals(this, obj)) return true;
            if (this is null || obj is null) return false;
            return Equals(obj as MElement);
        }

        public static MElement operator +(MElement a, MElement b)
        {
            if (a is null && b is null) return null;
            if (b is null) return a;
            if (a is null) return b;
            return a.AddTo(b);
        }

        public static MElement operator -(MElement a) => a?.Negate();
        public static MElement operator -(MElement a, MElement b) => a + -b;

        public static MElement operator *(MElement a, MElement b)
        {
            if (a is null || b is null) return null;
            return a.MultiplyTo(b);
        }

        public static MElement operator /(MElement a, MElement b) => a * b.Inverse();

        public static bool operator ==(MElement a, MElement b)
        {
            if (ReferenceEquals(a, b)) return true;
            if (a is null || b is null) return false;
            return a.Equals(b);
        }

        public static bool operator !=(MElement a, MElement b) => !(a == b);
    }
}
