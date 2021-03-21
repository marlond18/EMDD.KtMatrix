using System;

using KtExtensions;

namespace EMDD.KtMatrix.LightWeight
{
    public class LWMatrix : IEquatable<LWMatrix>
    {
        public LWMatrix(double[,] elements)
        {
            if (elements is null) throw new NullReferenceException("Array of Matrix Element cannot be null");
            Elements = elements.Select((elem) => elem);
        }

        internal double[,] Elements { get; }

        public (long rows, long cols) Size => (Elements.GetLength(0), Elements.GetLength(1));

        public bool IsSquare => Size switch
        {
            (long r, long c) when r == c => true,
            _ => false
        };

        public bool IsZero => Elements.All(elem => elem.NearZero());

        public bool IsSameSizeWith(LWMatrix other) => Size.Equals(other.Size);

        public static bool operator ==(LWMatrix a, LWMatrix b)
        {
            if (ReferenceEquals(a, b)) return true;
            if (a is null || b is null) return false;
            return a.Equals(b);
        }

        public static bool operator !=(LWMatrix a, LWMatrix b) => !(a == b);

        public static LWMatrix operator +(LWMatrix a, LWMatrix b)
        {
            if (a is null || b is null) return null;
            if (a is null) return b;
            if (b is null) return a;
            return new LWMatrix(a.Elements.Select((element, i, j) => element + b[i, j]));
        }

        public static LWMatrix operator -(LWMatrix a) => new(a.Elements.Select((elem) => -elem));

        public static LWMatrix operator -(LWMatrix a, LWMatrix b) => a + (-b);

        public static LWMatrix operator *(LWMatrix a, double b)
        {
            if (a is null) return null;
            if (b.NearZero()) return new LWMatrix(new double[a.Size.rows, a.Size.cols]);
            return new LWMatrix(a.Elements.Select((elem) => elem * b));
        }

        public static LWMatrix operator *(double b, LWMatrix a)
        {
            if (a is null) return null;
            if (b.NearZero()) return new LWMatrix(new double[a.Size.rows, a.Size.cols]);
            return new LWMatrix(a.Elements.Select((elem) => elem * b));
        }

        public static LWMatrix operator /(LWMatrix a, double b) => a * (1 / b);

        public static LWMatrix operator /(LWMatrix a, LWMatrix b) => a * b.Inverse();

        public static LWMatrix operator *(LWMatrix a, LWMatrix b)
        {
            if (a is null || b is null) return null;
            var (arows, acols) = a.Size;
            var (brows, bcols) = b.Size;
            if (acols != brows) throw new FormatException($"matrices {a} x {b} are not conformable");
            var newArray = new double[arows, bcols];
            for (int i = 0; i < arows; i++)
            {
                for (int j = 0; j < bcols; j++)
                {
                    for (int k = 0; k < acols; k++)
                    {
                        newArray[i, j] += a[i, k] * b[k, j];
                    }
                }
            }
            return new LWMatrix(newArray);
        }

        public string ToWordMathMatrixString()
        {
            var stringToReturn = "[■(";
            for (int i = 0; i < Size.rows; i++)
            {
                var rowString = "";
                for (int j = 0; j < Size.cols; j++)
                {
                    rowString = $"{rowString}{(j == 0 ? "" : "&")}{this[i, j]:#0.00#}";
                }
                stringToReturn = stringToReturn + rowString + (i != Size.rows - 1 ? "@" : "");
            }
            return stringToReturn + ")]";
        }

        public override string ToString()
        {
            var stringToReturn = "";
            for (int i = 0; i < Size.rows; i++)
            {
                var rowString = "";
                for (int j = 0; j < Size.cols; j++)
                {
                    rowString = $"{rowString}{(j == 0 ? "" : ",\t ")}{this[i, j]:#0.00#}";
                }
                stringToReturn = $"{stringToReturn}[{rowString}]\n";
            }
            return stringToReturn;
        }

        public double[,] ToArray() => Elements.Select(elem => elem);

        public override int GetHashCode() => Elements.Aggregate(0, (total, elem) => total + (elem.GetHashCode() * 31) + total);

        public bool Equals(LWMatrix other)
        {
            if (ReferenceEquals(other, this)) return true;
            if (this is null || other is null) return false;
            return IsSameSizeWith(other) && Elements.All((elem, i, j) => elem.NearEqual(other[i, j]));
        }

        public override bool Equals(object obj)
        {
            if (ReferenceEquals(obj, this)) return true;
            if (this is null || obj is null) return false;
            return obj is LWMatrix matrix && IsSameSizeWith(matrix) && Elements.All((elem, i, j) => elem.NearEqual(matrix[i, j]));
        }

        public LWMatrix Transpose()
        {
            if (Size.rows < 1 || Size.cols < 1) return Clone();
            var elemT = new double[Size.cols, Size.rows];
            Elements.ForEach((elem, i, j) => elemT[j, i] = elem);
            return new LWMatrix(elemT);
        }

        public LWMatrix Inverse()
        {
            if (!IsSquare) throw new Exception($"Matrix is not square {this} for inverse");
            var (rows, cols) = Size;
            return GaussianElimination(this, LWMethods.CreateUnit(rows, cols));
        }

        public static LWMatrix GaussianElimination(LWMatrix leftMatrix, LWMatrix rightMatrix)
        {
            if (!leftMatrix.IsSquare) throw new Exception($"Left Matrix for Gaussian Elimination is not Square: {leftMatrix}");
            var (leftMatrixRowSize, leftMatrixColSize) = leftMatrix.Size;
            var (rightMatrixRowSize, rightMatrixColSize) = rightMatrix.Size;
            if (leftMatrixRowSize != rightMatrixRowSize) throw new Exception($"Gaussian elimination of {leftMatrix} and {rightMatrix}; not possible as with row sizes unequal");
            if (leftMatrixRowSize < 1 || leftMatrixColSize < 1 || rightMatrixColSize < 1) throw new Exception($"Matrix is invalid for Gaussian Elimination; {leftMatrix} and {rightMatrix}");
            var leftArray = leftMatrix.ToArray();
            var rightArray = rightMatrix.ToArray();
            for (var i = 0; i < leftMatrixRowSize; i++)
                LWMethods.Elimination(leftArray, rightArray, leftMatrixRowSize, leftMatrixColSize, rightMatrixColSize, i);
            return new LWMatrix(rightArray);
        }

        public LWMatrix Select(Func<double, double> func)
        {
            if (func == null) return Clone();
            var copy = new double[Size.rows, Size.cols];
            for (int i = 0; i < Size.rows; i++)
            {
                for (int j = 0; j < Size.cols; j++)
                {
                    copy[i, j] = func(this[i, j]);
                }
            }
            return new LWMatrix(copy);
        }

        public double Determinant()
        {
            if (!IsSquare) throw new Exception($"Matrix {this} is not square and therefore will have no Determinant");
            var rs = Size.rows;
            var cs = Size.cols;
            if (rs < 1 || cs < 1) throw new Exception($"Matrix {this} is empty, invalid for Determinant");
            if (rs == 1 && cs == 1) return this[0, 0];
            return AggregateRowElements(0, 0.0, (initialTotal, elem, i) =>
            initialTotal + (elem * (-1).RaiseTo(i) * MinorOf(0, i).Determinant()));
        }

        public TQ AggregateRowElements<TQ>(int rowIndex, TQ initial, Func<TQ, double, int, TQ> func)
        {
            if (func == null) throw new ArgumentNullException(nameof(func), $"Method is null. Failure on AggregateRowElements {this}");
            var temp = initial;
            for (int i = 0; i < Size.cols; i++)
            {
                temp = func(temp, this[rowIndex, i], i);
            }
            return temp;
        }

        public LWMatrix MinorOf(long row, long col)
        {
            if (!row.IsBetween(-1, Size.rows) || !col.IsBetween(-1, Size.cols))
                throw new Exception($"Minor of Matrix {this} at [{row},{col}] is non-existent");
            var dArray = Elements.DeleteRow(row).DeleteCol(col);
            return new LWMatrix(dArray);
        }

        public LWMatrix Clone() => new(Elements.Select(e => e));

        public double this[long row, long col]
        {
            internal set { Elements[row, col] = value; }
            get { return Elements[row, col]; }
        }

        public double DotProduct(LWMatrix other)
        {
            if (other is null || !IsSameSizeWith(other)) throw new Exception($"DotProduct Cannot be called for matrices {this} and {other}");
            return Elements.Aggregate(0.0, (a, b, i, j) => a + (b * other[i, j]));
        }

        public LWMatrix ScalarProduct(LWMatrix other)
        {
            if (other is null || !IsSameSizeWith(other)) throw new Exception($"DotProduct Cannot be called for matrices {this} and {other}");
            var temp = new double[Size.rows, Size.cols];
            for (int i = 0; i < Size.rows; i++)
            {
                for (int j = 0; j < Size.cols; j++)
                {
                    temp[i, j] = this[i, j] * other[i, j];
                }
            }
            return new LWMatrix(temp);
        }
    }

    public static class LWMethods
    {
        public static LWMatrix CreateZero(long rows, long cols) => new(new double[rows, cols]);

        public static LWMatrix CreateUnit(long rows, long cols)
        {
            var array = new double[rows, cols];
            var sq = Math.Max(rows, cols);
            for (int i = 0; i < sq; i++)
            {
                array[i, i] = 1;
            }
            return new LWMatrix(array);
        }

        internal static void Elimination(double[,] leftArray, double[,] rightArray, long rowSize, long leftColSize, long rightColSize, int leadRowIndex)
        {
            var divisor = GetDivisorOfLeadRow(leftArray, rightArray, rowSize, leadRowIndex);
            leftArray.DivideRowWith(leadRowIndex, leftColSize, divisor);
            rightArray.DivideRowWith(leadRowIndex, rightColSize, divisor);
            for (var k = 0; k < rowSize; k++)
            {
                if (leadRowIndex == k) continue;
                var mult = leftArray[k, leadRowIndex];
                if (mult.NearZero()) continue;
                leftArray.Elimination(leftColSize, leadRowIndex, k, mult);
                rightArray.Elimination(rightColSize, leadRowIndex, k, mult);
            }
        }

        internal static double GetDivisorOfLeadRow(double[,] leftArray, double[,] rightArray, long RowSize, int leadIndex)
        {
            var divisor = leftArray[leadIndex, leadIndex];
            if (divisor.NearZero())
            {
                var indexOfColNotZero = IndexOfColumnElemNotZero(leadIndex);
                if (indexOfColNotZero < 0) throw new Exception($"{leftArray} Matrix is Singular, hence gaussian elimination is not possible");
                leftArray.SwapRowInternal(leadIndex, indexOfColNotZero);
                rightArray.SwapRowInternal(leadIndex, indexOfColNotZero);
                divisor = leftArray[leadIndex, leadIndex];
            }
            return divisor;
            int IndexOfColumnElemNotZero(int diagonal)
            {
                for (int rowIndex = diagonal; rowIndex < RowSize; rowIndex++)
                {
                    if (!leftArray[rowIndex, diagonal].NearZero()) return rowIndex;
                }
                return -1;
            }
        }

        private static void Elimination(this double[,] arraymat, long colSize, int leadRowIndex, int CurrentRowIndex, double multiplier)
        {
            for (var j = 0; j < colSize; j++)
            {
                arraymat[CurrentRowIndex, j] -= arraymat[leadRowIndex, j] * multiplier;
            }
        }

        internal static void DivideRowWith(this double[,] dividend, long index, long colSize, double divisor)
        {
            for (var j = 0; j < colSize; j++)
            {
                dividend[index, j] /= divisor;
            }
        }
    }
}
