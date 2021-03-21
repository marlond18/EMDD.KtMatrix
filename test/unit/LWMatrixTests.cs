using Microsoft.VisualStudio.TestTools.UnitTesting;
using LightWeightKtMatrix;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using EMDD.KtMatrix.LightWeight;

namespace LightWeightKtMatrix.Tests
{
    [TestClass()]
    public class LWMatrixTests
    {
        [TestMethod()]
        public void InverseTest()
        {
            var matrix = new LWMatrix(new[,] {
                {1.0 / 2, -2, -30, 20 },
                {15, 2.1, 4, 3},
                {3, 2, 1, 3   },
                { 4, -4, 3, 2 } });
            var determActual = matrix.Inverse();
            var determinantExpected = new LWMatrix(new[,]{
                { 0.0068357113280509, 0.097412395516784, -0.13109967926899, -0.017826187652207}, {-0.0092640030318555, 0.01172036747212, 0.15713714233579, -0.16064623439331}, { -0.027034045211142, -0.05670692764954, 0.18582747899808, 0.076659625088605}, {0.0083516390969001, -0.086323664615018, 0.29773242471243, 0.099370468884881 }});
            Assert.AreEqual(determinantExpected, determActual);
        }
    }
}