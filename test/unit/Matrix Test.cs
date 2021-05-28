using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using EMDD.KtMatrix;
using EMDD.KtPolynomials;
using System.Collections.Generic;
using System.Linq;
using EMDD.KtNumerics;
using EMDD.KtExpressions;
using KtExtensions;
using EMDD.KtExpressions.Limits;

namespace MatrixTest
{
    [TestClass]
    public partial class MatrixTest
    {
        [TestMethod]
        public void CastToZeroMatrixTest()
        {
            var matrix = new Matrix(new[,] {
                {1.0 / 2,    -2, -30,   20 },
                {     15,   2.1,   4,    3 },
                {      3,     2,   1,    3 },
                {      4,    -4,   3,    2 }
            });
            var zeros = matrix.ToZero();
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    Assert.AreEqual(zeros[i, j], 0);
                }
            }
        }

        [TestMethod]
        public void MatrixNotSquareExceptionTest()
        {
            var d = new double[,]
            {
                { 0.0,1.0,1},
                { 0,1,1}
            };
            var e = new double[,]
            {
                { 0.0,1.0,1},
                { 0,1,1}
            };
            var ex = Assert.ThrowsException<MatrixNotSquareExeption>(() =>
             {
                 var x = MatrixMethods.GaussianElimination(new Matrix(d), new Matrix(e));
             }
            );
            Console.WriteLine(ex.Message);
        }

        [TestMethod]
        public void CastToUnitMatrixTest()
        {
            var matrix = new Matrix(new[,] { {1.0 / 2, -2, -30, 20 },
                                                        {15, 2.1, 4, 3},
                                                        {3, 2, 1, 3   },
                                                        { 4, -4, 3, 2 } });
            var zeros = matrix.ToUnit();
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    if (i == j) Assert.AreEqual(zeros[i, j], 1);
                    else Assert.AreEqual(zeros[i, j], 0);
                }
            }
        }

        [TestMethod]
        public void DeterminantTest()
        {
            var matrix = new Matrix(new[,] { {1.0 / 2, -2, -30, 20 },
                                                        {15, 2.1, 4, 3},
                                                        {3, 2, 1, 3   },
                                                        { 4, -4, 3, 2 } });
            var determActual = matrix.Determinant();
            const double determinantExpected = -7124.35;
            Assert.AreEqual(determinantExpected, determActual);
        }

        [TestMethod]
        public void InverseTest()
        {
            var matrix = new Matrix(new[,] { {1.0 / 2, -2, -30, 20 },
                                                        {15, 2.1, 4, 3},
                                                        {3, 2, 1, 3   },
                                                        { 4, -4, 3, 2 } });
            var determActual = matrix.Inverse();
            var determinantExpected = new Matrix(new[,]{
                { 0.0068357113280509, 0.097412395516784, -0.13109967926899, -0.017826187652207}, {-0.0092640030318555, 0.01172036747212, 0.15713714233579, -0.16064623439331}, { -0.027034045211142, -0.05670692764954, 0.18582747899808, 0.076659625088605}, {0.0083516390969001, -0.086323664615018, 0.29773242471243, 0.099370468884881 }});
            Assert.AreEqual(determinantExpected, determActual);
        }
        [TestMethod]
        public void ComplexNumberTest()
        {
            var matrix = new EMDD.KtMatrix.Matrix(new EMDD.KtNumerics.Number[,] {
                { 3, 2,   1},
                { 3,  new KtComplex(3,1),   3},
                {-4, 1, 0.4} });
            var eq = EMDD.KtPolynomials.KtPolynomial.Create(3, -2);
            var lim = EMDD.KtExpressions.Limits.Limit.Create(3, 4);
            var exp = new EMDD.KtExpressions.Expression((eq, lim));
            var comp = new EMDD.KtNumerics.KtComplex(3, 1);

            var matrix2 = new EMDD.KtMatrix.Matrix(new EMDD.KtExpressions.Expression[,] {
        { 3,       2,    1},
        { exp,  comp,    3},
        {-4,       1,  0.4}
            });
        }

        [TestMethod]
        public void MultiplicationOfExpressionAndNumber()
        {
            var mat1 = new Matrix(new[,] {
                {1.0 / 2,       -2 ,        -30,     20},
                {15     ,       2.1,        4  ,     3 },
                {3      ,       2  ,        1  ,     3 },
                { 4     ,       -4 ,        3  ,     2 }});
            var expression1 = new Expression(((KtPolynomial.Create(3, 2), KtPolynomial.Create(1, 2)), Limit.Create(3, 12)));
            var expression2 = new Expression(((KtPolynomial.Create(2), KtPolynomial.Create(1, 2)), Limit.Create(3, 12)));

            var mat2 = new Matrix(new[,] {
                {expression1},{ expression2},{expression1},{ expression2} });
            _ = mat1 * mat2;
        }

        [TestMethod]
        public void MultiplicationTest()
        {
            var Beams = new List<double>() { 1000, 500, 500 };
            var BeamSpanDistances = Beams;
            var L = BeamSpanDistances.Sum();
            var beamCount = Beams.Count;
            var Coeffs = new Number[beamCount - 1, beamCount - 1];
            var Expressions = new Expression[beamCount - 1, 1];
            var xVar = Monomial.Create(1, 1);
            var CommulativeLengths = BeamSpanDistances.AggregateSelect(0.0, (previousCollective, current) => previousCollective + current).ToArray();
            for (int i = 0; i < beamCount - 1; i++)
            {
                var a = xVar;
                var b = L - xVar;
                var supportLocation = CommulativeLengths[i + 1];
                var x = Monomial.Create(0, supportLocation);
                Expressions[i, 0] = RightSideExpression(a, x, 0, supportLocation, L) + LeftSideExpression(a, x, supportLocation, L, L);
                for (int j = 0; j < beamCount - 1; j++)
                {
                    var aa = CommulativeLengths[j];
                    var bb = L - aa;
                    var xx = supportLocation;
                    Coeffs[i, j] = (j < i) ? RightSideExpression(aa, bb).Evaluate(xx) : LeftSideExpression(aa, bb).Evaluate(xx);
                }
            }
            var mat1 = new Matrix(Coeffs);
            var mat2 = new Matrix(Expressions);
            var inverse = mat1.Inverse();
            var mat3 = (Matrix)(inverse * mat2);
        }

        private static Expression LeftSideExpression(double a, double b)
        {
            var L = a + b;
            var kteq = KtPolynomial.Create(-a.RaiseTo(3) / 6, (a.RaiseTo(3) / (6 * L)) + ((L * a) / 3), -a / 2, a / (6 * L));
            return new Expression(kteq);
        }

        private static Expression LeftSideExpression(KtPolynomial a, KtPolynomial x, double limit1, double limit2, double L)
        {
            var kteq = (a - L) * x * ((x * x) - (2 * a * L) + (a * a)) / (6 * L);
            return new Expression((kteq, Limit.Create(limit1, limit2)));
        }

        private static Expression RightSideExpression(double a, double b)
        {
            var L = a + b;
            var kteq = KtPolynomial.Create(0, ((a * b * b) / (3 * L)) + ((a * a * b) / (6 * L)), 0, -b / (6 * L));
            return new Expression(kteq);
        }

        private static Expression RightSideExpression(KtPolynomial a, KtPolynomial x, double limit1, double limit2, double L)
        {
            var kteq = a * (x - L) * ((x * x) - (2 * L * x) + (a * a)) / (6 * L);
            return new Expression((kteq, Limit.Create(limit1, limit2)));
        }

        private const int interval = 2000;

#pragma warning disable CA1822 // Mark members as static
#pragma warning disable IDE0051 // Remove unused private members
        private void RunGauElim(Matrix mat1, Matrix mat2)
#pragma warning restore IDE0051 // Remove unused private members
#pragma warning restore CA1822 // Mark members as static
        {
            for (int i = 0; i < interval; i++)
            {
                _ = MatrixMethods.GaussianElimination(mat1, mat2);
            }
        }

        [TestMethod]
        public void GaussianEliminationTest()
        {
            var mat1 = new Matrix(new[,]{
                { 1936.0 / 45, 692.0 / 15},
                { 692.0 / 15, 324.0 / 5 } });
            var polyu = KtPolynomial.Create(1);
            var poly1 = KtPolynomial.Create(0, 572.0 / 45, 0, -11.0 / 90);
            var poly2 = KtPolynomial.Create(-32.0 / 3, 932.0 / 45, -2, 2.0 / 45);
            var expression1 = new Expression(
                ((poly1, polyu), Limit.Create(0, 4)),
                ((poly2, polyu), Limit.Create(4, 9)),
                ((poly2, polyu), Limit.Create(9, 15)));
            var poly3 = KtPolynomial.Create(0, 63.0 / 5, 0, -1 / 15.0);
            var poly4 = KtPolynomial.Create(-243 / 2.0, 531.0 / 10, -9.0 / 2, 1.0 / 10);
            var expression2 = new Expression(
                ((poly3, polyu), Limit.Create(0, 4)),
                ((poly3, polyu), Limit.Create(4, 9)),
                ((poly4, polyu), Limit.Create(9, 15)));
            var mat2 = new Matrix(new[,] { { expression1 }, { expression2 } });

            var actual = MatrixMethods.GaussianElimination(mat1, mat2);
            //GetEllapsedTime(mat1, mat2);
            var expPoly1a = KtPolynomial.Create(0, 2727.0 / 7420, 0, -109.0 / 14840);
            var expPoly1b = KtPolynomial.Create(-1944.0 / 1855, 8559.0 / 7420, -729.0 / 3710, 67.0 / 7420);
            var expPoly1c = KtPolynomial.Create(3159.0 / 424, -24921.0 / 14840, 351.0 / 2968, -39.0 / 14840);
            var expPoly2a = KtPolynomial.Create(0, -374.0 / 5565, 0, 187.0 / 44520);
            var expPoly2b = KtPolynomial.Create(1384.0 / 1855, -3488.0 / 5565, 519.0 / 3710, -83.0 / 11130);
            var expPoly2c = KtPolynomial.Create(-761.0 / 106, 22427.0 / 11130, -57.0 / 371, 19.0 / 5565);
            var expected1 = new Expression(
                ((expPoly1a, polyu), Limit.Create(0, 4)),
                ((expPoly1b, polyu), Limit.Create(4, 9)),
                ((expPoly1c, polyu), Limit.Create(9, 15)));
            var expected2 = new Expression(
                ((expPoly2a, polyu), Limit.Create(0, 4)),
                ((expPoly2b, polyu), Limit.Create(4, 9)),
                ((expPoly2c, polyu), Limit.Create(9, 15)));
            var expected = new Matrix(new[,] { { expected1 }, { expected2 } });
            Assert.AreEqual(expected1, actual[0, 0]);
            Assert.AreEqual(expected2, actual[1, 0]);
            Assert.AreEqual(expected, actual);
        }

#pragma warning disable IDE0051 // Remove unused private members
        private static void GetEllapsedTime(Matrix mat1, Matrix mat2)
#pragma warning restore IDE0051 // Remove unused private members
        {
            Console.WriteLine(DiagnosticsExtensions.ElapsedTime(() => { for (int i = 0; i < interval; i++) { var actual2 = MatrixMethods.GaussianElimination(mat1, mat2); } }));
            Console.WriteLine(DiagnosticsExtensions.ElapsedTime(() => { for (int i = 0; i < interval; i++) { var actual2 = MatrixMethods.GaussianElimination2(mat1, mat2); } }));
            Console.WriteLine(DiagnosticsExtensions.ElapsedTime(() => { for (int i = 0; i < interval; i++) { var actual2 = MatrixMethods.GaussianElimination2(mat1, mat2); } }));
            Console.WriteLine(DiagnosticsExtensions.ElapsedTime(() => { for (int i = 0; i < interval; i++) { var actual1 = MatrixMethods.GaussianElimination(mat1, mat2); } }));
            Console.WriteLine(DiagnosticsExtensions.ElapsedTime(() => { for (int i = 0; i < interval; i++) { var actual1 = MatrixMethods.GaussianElimination(mat1, mat2); } }));
            Console.WriteLine(DiagnosticsExtensions.ElapsedTime(() => { for (int i = 0; i < interval; i++) { var actual2 = MatrixMethods.GaussianElimination2(mat1, mat2); } }));
            Console.WriteLine(DiagnosticsExtensions.ElapsedTime(() => { for (int i = 0; i < interval; i++) { var actual2 = MatrixMethods.GaussianElimination2(mat1, mat2); } }));
            Console.WriteLine(DiagnosticsExtensions.ElapsedTime(() => { for (int i = 0; i < interval; i++) { var actual1 = MatrixMethods.GaussianElimination(mat1, mat2); } }));
            Console.WriteLine(DiagnosticsExtensions.ElapsedTime(() => { for (int i = 0; i < interval; i++) { var actual1 = MatrixMethods.GaussianElimination(mat1, mat2); } }));
            Console.WriteLine(DiagnosticsExtensions.ElapsedTime(() => { for (int i = 0; i < interval; i++) { var actual2 = MatrixMethods.GaussianElimination2(mat1, mat2); } }));
            Console.WriteLine(DiagnosticsExtensions.ElapsedTime(() => { for (int i = 0; i < interval; i++) { var actual2 = MatrixMethods.GaussianElimination2(mat1, mat2); } }));
            Console.WriteLine(DiagnosticsExtensions.ElapsedTime(() => { for (int i = 0; i < interval; i++) { var actual1 = MatrixMethods.GaussianElimination(mat1, mat2); } }));
        }
    }
}