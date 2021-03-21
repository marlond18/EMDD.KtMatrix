using System;
using EMDD.KtNumerics;
using EMDD.KtExpressions;
using System.Collections.Generic;
using System.Linq;
using KtExtensions;

namespace EMDD.KtMatrix
{
    public class Matrix : MElement
    {
        public Matrix(double[,] elements)
        {
            if (elements is null) throw new NullReferenceException("Array of Matrix Element cannot be null");
            Elements = elements.Select((elem) => (MElement)elem);
        }

        public Matrix(Number[,] elements)
        {
            if (elements is null) throw new NullReferenceException("Array of Matrix Element cannot be null");
            Elements = elements.Select((elem) => (MElement)elem);
        }

        public Matrix(Expression[,] elements)
        {
            if (elements is null) throw new NullReferenceException("Array of Matrix Element cannot be null");
            Elements = elements.Select((elem) => (MElement)elem);
        }

        public Matrix(MElement[,] elements)
        {
            if (elements is null) throw new NullReferenceException("Array of Matrix Element cannot be null");
            Elements = elements.Select((elem) => elem.Clone());
        }

        public Matrix Col(int index)
        {
            if (index < 0 || index > Size.cols - 1) return null;
            var col = new MElement[Size.rows, 1];
            for (int i = 0; i < Size.rows; i++)
            {
                col[i, 0] = this[i, index];
            }
            return new Matrix(col);
        }

        public Matrix Row(int index)
        {
            if (index < 0 || index > Size.rows - 1) return null;
            var row = new MElement[1, Size.cols];
            for (int i = 0; i < Size.cols; i++)
            {
                row[0, i] = this[index, i];
            }
            return new Matrix(row);
        }

        internal MElement[,] Elements { get; }

        public (long rows, long cols) Size => (Elements.GetLength(0), Elements.GetLength(1));

        public bool IsSquare => Size.rows == Size.cols;

        public override bool IsZero => Elements.All(elem => elem.IsZero);

        public bool IsSameSizeWith(Matrix other) => Size.Equals(other.Size);

        public override MElement AddTo(MElement other)
        {
            if (other is Matrix mat)
            {
                if (!IsSameSizeWith(mat)) throw new Exception("Matrix cannot be added to this number since they are of different size");
                return new Matrix(Elements.Select((element, i, j) => element + mat[i, j]));
            }
            return this + new Matrix(new[,] { { other } });
        }

        public override MElement Negate() => new Matrix(Elements.Select((elem) => elem.Negate()));

        public Matrix MultiplyTo(Matrix matrix)
        {
            if (Size.cols != matrix.Size.rows) throw new FormatException($"matrices {this} x {matrix} are not conformable");
            var newArray = new MElement[Size.rows, matrix.Size.cols];
            for (int i = 0; i < Size.rows; i++)
            {
                for (int j = 0; j < matrix.Size.cols; j++)
                {
                    newArray[i, j] = SelectRow(i).Zip(matrix.SelectCol(j), (e1, e2) => e1 * e2).Sum();
                }
            }
            return new Matrix(newArray);
        }

        public override MElement MultiplyTo(MElement other) => other switch
        {
            MNumeric numeric => new Matrix(Elements.Select((elem) => elem * numeric._value)),
            MExpression expression => new Matrix(Elements.Select((elem) => elem * expression._value)),
            Matrix matrix => MultiplyTo(matrix),
            _ => throw new Exception($"{nameof(MElement)} - ({other.GetType()}) type is invalid to be multiplied to a Matrix"),
        };

        public string ToWordMathMatrixString() => $"[■({SelectRows(enumerable => enumerable.BuildString("&").BuildString("@"))})]";

        public string ToWordMathMatrixStringPieceWise() => $"[■({SelectRows(element => element is MExpression exp ? exp.ToStringPieceWise() : element.ToString(), enumerable => enumerable.BuildString("&").BuildString("@"))})]";

        public override string ToString() => SelectRows(e => e.ToString(), s => s.BuildString(", ")).BuildString(",\n", ("{", "}"));

        public MElement[,] ToArray() => Elements.Select((elem) => elem.Clone());

        public override int GetHashCode() => Elements.Aggregate(0, (total, elem) => total + (elem.GetHashCode() * 31) + total);

        public override bool Equals(MElement other) => other is Matrix matrix && IsSameSizeWith(matrix) && Elements.All((elem, i, j) => elem == matrix[i, j]);

        public Matrix Transpose()
        {
            if (Size.rows < 1 || Size.cols < 1) return Clone() as Matrix;
            var elemT = new MElement[Size.rows, Size.cols];
            Elements.ForEach((elem, i, j) => elemT[j, i] = elem);
            return new Matrix(elemT);
        }

        public override MElement Inverse()
        {
            if (!IsSquare) throw new Exception($"Matrix is not square {this} for inverse");
            return MatrixMethods.GaussianElimination(this, ToUnit());
        }

        public Matrix ToUnit() => new(Elements.Select((MElement elem, int i, int j) => (elem is Matrix matrix)
            ? i == j ? matrix.ToUnit() : matrix.ToZero()
            : i == j ? (MElement)1 : 0));

        public Matrix ToZero() => new(Elements.Select(elem => (elem is Matrix matrix) ? matrix.ToZero() : (MElement)0));

        public MElement Determinant()
        {
            if (!IsSquare) throw new Exception($"Matrix {this} is not square and therefore will have no Determinant");
            var (rs, cs) = Size;
            if (rs < 1 || cs < 1) throw new Exception($"Matrix {this} is empty, invalid for Determinant");
            if (rs == 1 && cs == 1) return this[0, 0];
            return AggregateRowElements(0, (MElement)0, (initialTotal, elem, i) => initialTotal + (elem * (-1).RaiseTo(i) * MinorOf(0, i).Determinant()));
        }

        public TQ AggregateRowElements<TQ>(int rowIndex, TQ initial, Func<TQ, MElement, int, TQ> func)
        {
            if (func == null) throw new ArgumentNullException(nameof(func), $"Method is null. Failure on AggregateRowElements {this}");
            var temp = initial;
            for (int i = 0; i < Size.cols; i++)
            {
                temp = func(temp, this[rowIndex, i], i);
            }
            return temp;
        }

        public Matrix MinorOf(long row, long col)
        {
            if (!row.IsBetween(-1, Size.rows) || !col.IsBetween(-1, Size.cols)) throw new Exception($"Minor of Matrix {this} at [{row},{col}] is non-existent");
            return new Matrix(Elements.DeleteRow(row).DeleteCol(col));
        }

        public override MElement Clone() => new Matrix(Elements.Select(e => e.Clone()));

        public MElement this[long row, long col]
        {
            get => Elements[row, col];
            internal set { Elements[row, col] = value; }
        }

        public MElement DotProduct(Matrix other)
        {
            if (other == null || !IsSameSizeWith(other)) throw new Exception($"DotProduct Cannot be called for matrices {this} and {other}");
            return Elements.Aggregate((MElement)0, (a, b, i, j) => a + (b * other[i, j]));
        }

        public override double ToDouble() => 0;

        public IEnumerable<T> SelectCol<T>(int index, Func<MElement, T> func) => Elements.SelectCol(index, func);

        public IEnumerable<MElement> SelectCol(int index) => Elements.SelectCol(index);

        public IEnumerable<T> SelectRow<T>(int index, Func<MElement, T> func) => Elements.SelectRow(index, func);

        public IEnumerable<MElement> SelectRow(int index) => Elements.SelectRow(index);

        public IEnumerable<T> SelectRows<T>(Func<IEnumerable<MElement>, T> func) => Elements.SelectRows(func);

        public IEnumerable<T2> SelectRows<T1, T2>(Func<MElement, T1> func, Func<IEnumerable<T1>, T2> func2) => Elements.SelectRows(func, func2);

        public IEnumerable<T> SelectCols<T>(Func<IEnumerable<MElement>, T> func) => Elements.SelectCols(func);

        public IEnumerable<T2> SelectCols<T1, T2>(Func<MElement, T1> func, Func<IEnumerable<T1>, T2> func2) => Elements.SelectCols(func, func2);
    }
}