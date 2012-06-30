
#include "../cute/cute.h"
#include "../cute/ide_listener.h"
#include "../cute/cute_suite.h"

#include "../partyql/graph/Graph.h"
#include "../partyql/basicmath/scalar/Integer.h"
#include "../partyql/basicmath/Vektor.h"
#include "../partyql/basicmath/Matrix.h"

using namespace AlgoTrans;

namespace Test_Matrix {

	CInteger aiv4_4_1[] = {CInteger(4), CInteger(4), CInteger(1)};
	CInteger aiv2__3__1[] = {CInteger(2), CInteger(-3), CInteger(-1)};
	CInteger aiv4_3_2[] = {CInteger(4), CInteger(3), CInteger(2)};

	CInteger aiv6_0_1[] = {CInteger(6), CInteger(0), CInteger(1)};

	CInteger aiv0_1__1[] = {CInteger(0), CInteger(1), CInteger(-1)};
	CInteger aiv0_0_13[] = {CInteger(0), CInteger(0), CInteger(13)};

	CInteger aiv2_0_9[] = {CInteger(2), CInteger(0), CInteger(9)};
	CInteger aiv0_1_12[] = {CInteger(0), CInteger(1), CInteger(12)};

	CInteger aiv3_8__18[] = {CInteger(3), CInteger(8), CInteger(-18)};

	CInteger aiv1_0_0_1[]  = {CInteger(1), CInteger(0), CInteger(0), CInteger(1)};
	CInteger aiv0_1_0_0[]  = {CInteger(0), CInteger(1), CInteger(0), CInteger(0)};
	CInteger aiv0_0_1_0[]  = {CInteger(0), CInteger(0), CInteger(1), CInteger(0)};
	CInteger aiv0_1__1_0[] = {CInteger(0), CInteger(1), CInteger(-1), CInteger(0)};
	CInteger aiv1_0_0__1[]  = {CInteger(1), CInteger(0), CInteger(0), CInteger(-1)};
	CInteger aiv0_1_1_0[]  = {CInteger(0), CInteger(1), CInteger(1), CInteger(0)};

	CInteger aiv0__1[]  = {CInteger(0), CInteger(-1)};
	CInteger aiv1_0[]  = {CInteger(1), CInteger(0)};
	CInteger aiv_1_0[]  = {CInteger(-1), CInteger(0)};

	CVector<CInteger> iv0__1, iv1_0, iv_1_0;

	CVector<CInteger> iv4_4_1, iv4_3_2, iv6_0_1;
	CVector<CInteger> iv2__3__1, iv0_1__1, iv0_0_13;
	CVector<CInteger> iv2_0_9, iv0_1_12;

	CVector<CInteger> iv3_8__18;

	CVector<CInteger> iv1_0_0_1, iv0_1_0_0, iv0_0_1_0, iv0_1__1_0, iv1_0_0__1, iv0_1_1_0;

	void test_Matrix_initializeIntegerVectors() {
		iv4_4_1 = CVector<CInteger>(3, aiv4_4_1);
		iv2__3__1 = CVector<CInteger>(3, aiv2__3__1);
		iv6_0_1 = CVector<CInteger>(3, aiv6_0_1);
		iv4_3_2 = CVector<CInteger>(3, aiv4_3_2);

		iv2__3__1 = CVector<CInteger>(3, aiv2__3__1);
		iv0_1__1 = CVector<CInteger>(3, aiv0_1__1);
		iv0_0_13 = CVector<CInteger>(3, aiv0_0_13);

		iv2_0_9 = CVector<CInteger>(3, aiv2_0_9);
		iv0_1_12 = CVector<CInteger>(3, aiv0_1_12);
		iv3_8__18 = CVector<CInteger>(3, aiv3_8__18);

		iv1_0_0_1 = CVector<CInteger>(4, aiv1_0_0_1);
		iv0_1_0_0 = CVector<CInteger>(4, aiv0_1_0_0);
		iv0_0_1_0 = CVector<CInteger>(4, aiv0_0_1_0);
		iv0_1__1_0 = CVector<CInteger>(4, aiv0_1__1_0);
		iv1_0_0__1 = CVector<CInteger>(4, aiv1_0_0__1);
		iv0_1_1_0 = CVector<CInteger>(4, aiv0_1_1_0);

		iv0__1 = CVector<CInteger>(2, aiv0__1);
		iv1_0 = CVector<CInteger>(2, aiv1_0);
		iv_1_0 = CVector<CInteger>(2, aiv_1_0);
	}

	void CMatrix_HermiteNormalFormReduction() {
		test_Matrix_initializeIntegerVectors();

		CMatrix<CInteger> a = CMatrix<CInteger>(3);
		a.addRow(iv4_4_1);
		a.addRow(iv2__3__1);
		a.addRow(iv4_3_2);
		a.reduceToHermiteNormalForm();

		CMatrix<CInteger> exp = CMatrix<CInteger>(3);
		exp.addRow(iv2_0_9);
		exp.addRow(iv0_1_12);
		exp.addRow(iv0_0_13);

		ASSERT(a == exp);

		// This one should have the same solution
		CMatrix<CInteger> c = CMatrix<CInteger>(3);
		c.addRow(iv4_4_1);
		c.addRow(iv6_0_1);
		c.addRow(iv4_3_2);
		c.reduceToHermiteNormalForm();

		ASSERT(c == exp); // XXX: Err... we should test whether matrix equality works too :)
	}

}

cute::suite* Test_Matrix_runSuite(){
	cute::suite &s = *(new cute::suite("Matrix"));

	s.push_back(CUTE(Test_Matrix::CMatrix_HermiteNormalFormReduction));

	return &s;
}


