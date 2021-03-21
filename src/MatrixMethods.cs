using System;
using System.Collections.Generic;
using System.Linq;
using System.Collections;
using System.Runtime.CompilerServices;
using KtExtensions;

namespace EMDD.KtMatrix
{
    public class MatrixNotSquareExeption : Exception
    {
        private static string CreateMessage(Matrix matrix, string methodName) => $"{matrix} used in {methodName} is not square, operation is terminated";

        public MatrixNotSquareExeption(Matrix matrix, [CallerMemberName] string methodName = null) : base(CreateMessage(matrix, methodName))
        {
        }

        protected MatrixNotSquareExeption(string message) : base(message)
        {
        }

        protected MatrixNotSquareExeption(string message, Exception innerException) : base(message, innerException)
        {
        }

        protected MatrixNotSquareExeption(System.Runtime.Serialization.SerializationInfo info, System.Runtime.Serialization.StreamingContext context) : base(info, context)
        {
        }

        protected MatrixNotSquareExeption() : base()
        {
        }
    }

    public class MatrixRow : IEquatable<MatrixRow>, IEnumerable<MElement>
    {
        private readonly MElement[] _elements;

        public MatrixRow(IEnumerable<MElement> elements)
        {
            _elements = elements.ToArray();
            Size = _elements.Length;
        }

        public int Size { get; }

        public bool Equals(MatrixRow other)
        {
            if (ReferenceEquals(this, other)) return true;
            if (this is null || other is null) return false;
            return Size == other.Size && _elements.SequenceEqual(other._elements);
        }

        public override int GetHashCode() => _elements.GetHashCodeOfEnumerable();

        public override bool Equals(object obj) => Equals(obj as MatrixRow);

        public IEnumerator<MElement> GetEnumerator() => (IEnumerator<MElement>)_elements.GetEnumerator();

        IEnumerator IEnumerable.GetEnumerator() => _elements.GetEnumerator();

        public static MatrixRow operator +(MatrixRow a, MatrixRow b) => new(a._elements.Zip(b._elements, (e1, e2) => e1 + e2));

        public static MatrixRow operator -(MatrixRow a, MatrixRow b) => new(a._elements.Zip(b._elements, (e1, e2) => e1 - e2));

        public static MatrixRow operator -(MatrixRow a) => new(a._elements.Select(e => -e));

        public static MatrixRow operator /(MatrixRow a, MElement b) => new(a._elements.Select(e => e / b));

        public static MatrixRow operator *(MatrixRow a, MElement b) => new(a._elements.Select(e => e * b));

        public static bool operator ==(MatrixRow row1, MatrixRow row2) => EqualityComparer<MatrixRow>.Default.Equals(row1, row2);

        public static bool operator !=(MatrixRow row1, MatrixRow row2) => !(row1 == row2);

        public MElement this[int i] => _elements[i];
    }

    public static class MatrixMethods
    {
        internal static IEnumerable<MatrixRow> ToRows(this Matrix mat) => mat.SelectRows(ienum => new MatrixRow(ienum));

        internal static Matrix CreateMatrix(this IEnumerable<MatrixRow> rows)
        {
            var (colCount, rowCount) = (rows.Min(row => row.Size), rows.Count());
            var mat = new MElement[rowCount, colCount];
            var list = rows.ToList();
            for (int i = 0; i < rowCount; i++)
            {
                for (int j = 0; j < colCount; j++)
                {
                    mat[i, j] = list[i][j];
                }
            }
            return new Matrix(mat);
        }

        public static MElement Sum(this IEnumerable<MElement> enumerable) => enumerable.Aggregate((MElement)0, (total, e) => total + e);

        public static Matrix CreateZero(int row, int col) => new(new double[row, col]);

        public static Matrix GaussianElimination(Matrix leftMatrix, Matrix rightMatrix)
        {
            if (!leftMatrix.IsSquare) throw new MatrixNotSquareExeption(leftMatrix);
            return GaussianInner(leftMatrix.ToRows().ToList(), rightMatrix.ToRows().ToList(), GetRowSize(leftMatrix, rightMatrix));
        }

        private static Matrix GaussianInner(List<MatrixRow> leftRows, List<MatrixRow> rightRows, long rowSize)
        {
            for (var i = 0; i < rowSize; i++)
            {
                NormalizeLeadRow(leftRows, rightRows, i, GetDivisor(leftRows, rightRows, i));
                EliminateColumn(leftRows, rightRows, i, rowSize);
            }
            return rightRows.CreateMatrix();
        }

        private static long GetRowSize(Matrix leftMatrix, Matrix rightMatrix)
        {
            var (leftMatrixRowSize, leftMatrixColSize) = leftMatrix.Size;
            var (rightMatrixRowSize, rightMatrixColSize) = rightMatrix.Size;
            if (leftMatrixRowSize != rightMatrixRowSize) throw new Exception($"Gaussian elimination of {leftMatrix} and {rightMatrix}; not possible as with row sizes unequal");
            if (leftMatrixRowSize < 1 || leftMatrixColSize < 1 || rightMatrixColSize < 1) throw new Exception($"Matrix is invalid for Gaussian Elimination; {leftMatrix} and {rightMatrix}");
            return leftMatrixRowSize;
        }

        private static void NormalizeLeadRow(List<MatrixRow> leftRows, List<MatrixRow> rightRows, int lead, MElement divisor)
        {
            leftRows[lead] /= divisor;
            rightRows[lead] /= divisor;
        }

        private static MElement GetDivisor(List<MatrixRow> leftRows, List<MatrixRow> rightRows, int lead)
        {
            if (!leftRows[lead][lead].IsZero) return leftRows[lead][lead];
            SwapNoneZeroLead(leftRows, rightRows, lead);
            return leftRows[lead][lead];
        }

        private static void SwapNoneZeroLead(List<MatrixRow> leftRows, List<MatrixRow> rightRows, int lead)
        {
            int indexOfColNotZero = leftRows.FindIndex(row => !row[lead].IsZero);
            if (indexOfColNotZero < 0) throw new Exception("Left Matrix is Singular, hence gaussian elimination is not possible");
            leftRows.Swap(lead, indexOfColNotZero);
            rightRows.Swap(lead, indexOfColNotZero);
        }

        private static void EliminateColumn(List<MatrixRow> leftRows, List<MatrixRow> rightRows, int column, long rowSize)
        {
            for (var k = 0; k < rowSize; k++)
            {
                if (column == k) continue;
                var mult = leftRows[k][column];
                if (mult.IsZero) continue;
                leftRows[k] -= (leftRows[column] * mult);
                rightRows[k] -= (rightRows[column] * mult);
            }
        }

        public static Matrix GaussianElimination2(Matrix leftMatrix, Matrix rightMatrix)
        {
            if (!leftMatrix.IsSquare) throw new Exception($"Left Matrix for Gaussian Elimination is not Square: {leftMatrix}");
            var (leftMatrixRowSize, leftMatrixColSize) = leftMatrix.Size;
            var (rightMatrixRowSize, rightMatrixColSize) = rightMatrix.Size;
            if (leftMatrixRowSize != rightMatrixRowSize) throw new Exception($"Gaussian elimination of {leftMatrix} and {rightMatrix}; not possible as with row sizes unequal");
            if (leftMatrixRowSize < 1 || leftMatrixColSize < 1 || rightMatrixColSize < 1) throw new Exception($"Matrix is invalid for Gaussian Elimination; {leftMatrix} and {rightMatrix}");
            var leftArray = leftMatrix.ToArray();
            var rightArray = rightMatrix.ToArray();
            for (var i = 0; i < leftMatrixRowSize; i++)
            {
                var divisor = leftArray[i, i];
                if (divisor.IsZero)
                {
                    int indexOfColNotZero = leftArray.NoneZeroOfColumn(leftMatrixRowSize, i);
                    if (indexOfColNotZero < 0) throw new Exception("Left Matrix is Singular, hence gaussian elimination is not possible");
                    leftArray.SwapRowInternal(i, indexOfColNotZero);
                    rightArray.SwapRowInternal(i, indexOfColNotZero);
                    divisor = leftArray[i, i];
                }
                leftArray.DivideRowWith(i, leftMatrixColSize, divisor);
                rightArray.DivideRowWith(i, rightMatrixColSize, divisor);
                for (var k = 0; k < leftMatrixRowSize; k++)
                {
                    if (i == k) continue;
                    var mult = leftArray[k, i];
                    if (mult.IsZero) continue;
                    leftArray.Elimination(leftMatrixColSize, i, k, mult);
                    rightArray.Elimination(rightMatrixColSize, i, k, mult);
                }
            }
            return new Matrix(rightArray);
        }

        internal static void DivideRowWith(this MElement[,] dividend, long index, long colSize, MElement divisor)
        {
            for (var j = 0; j < colSize; j++)
            {
                dividend[index, j] /= divisor;
            }
        }

        internal static int NoneZeroOfColumn(this MElement[,] leftArray, long rowSize, int columnIndex)
        {
            for (int rowIndex = columnIndex; rowIndex < rowSize; rowIndex++)
            {
                if (!leftArray[rowIndex, columnIndex].IsZero) return rowIndex;
            }
            return -1;
        }

        internal static void Elimination(this MElement[,] arraymat, long colSize, int leadRowIndex, int CurrentRowIndex, MElement multiplier)
        {
            for (var j = 0; j < colSize; j++)
            {
                arraymat[CurrentRowIndex, j] -= arraymat[leadRowIndex, j] * multiplier;
            }
        }
    }
}
