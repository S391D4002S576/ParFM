#include "../cute/cute.h"
#include "../cute/ide_listener.h"
#include "../cute/cute_suite.h"

#include "Test_Funcs.h"

#include "../partyql/graph/Graph.h"
#include "../partyql/basicmath/scalar/Integer.h"
#include "../partyql/basicmath/Vektor.h"
#include "../partyql/basicmath/Module.h"
#include "../partyql/basicmath/SetRelation.h"
#include "../partyql/basicmath/Hadron.h"
#include "../partyql/basicmath/VectorDescriptor.h"

using namespace AlgoTrans;

namespace Test_ModuleRelation {
	CInteger aiv1_0_0_1[] = {CInteger(1), CInteger(0), CInteger(0), CInteger(1)};
	CInteger aiv0_1_0_0[] = {CInteger(0), CInteger(1), CInteger(0), CInteger(0)};
	CInteger aiv0_0_1_0[] = {CInteger(0), CInteger(0), CInteger(1), CInteger(0)};
	CVector<CInteger> iv1_0_0_1, iv0_1_0_0, iv0_0_1_0;

	CInteger aiv1_0_0_0_0_1[] = {CInteger(1), CInteger(0), CInteger(0), CInteger(0), CInteger(0), CInteger(1)};
	CInteger aiv0_1_0_0_0_0[] = {CInteger(0), CInteger(1), CInteger(0), CInteger(0), CInteger(0), CInteger(0)};
	CInteger aiv0_0_1_0_0_0[] = {CInteger(0), CInteger(0), CInteger(1), CInteger(0), CInteger(0), CInteger(0)};
	CVector<CInteger> iv1_0_0_0_0_1, iv0_1_0_0_0_0, iv0_0_1_0_0_0;

	CInteger aiv1_0_0__1[] = {CInteger(1), CInteger(0), CInteger(0),  CInteger(-1)};
	CInteger aiv0_1__1_0[] = {CInteger(0), CInteger(1), CInteger(-1), CInteger(0)};
	CVector<CInteger> iv1_0_0__1, iv0_1__1_0;

	CInteger aiv0_0_0__1[] = {CInteger(0), CInteger(0), CInteger(0),  CInteger(-1)};
	CInteger aiv1_0__1_0[] = {CInteger(1), CInteger(0), CInteger(-1), CInteger(0)};
	CVector<CInteger> iv0_0_0__1, iv1_0__1_0;

	CInteger aiv0_1_1_1[] = {CInteger(0), CInteger(1), CInteger(1), CInteger(1)};
	CInteger aiv0_1__1__1[] = {CInteger(0), CInteger(1), CInteger(-1), CInteger(-1)};
	CVector<CInteger> iv0_1_1_1, iv0_1__1__1;

	CInteger aiv1_0_0_1_11_7[] = {CInteger(1), CInteger(0), CInteger(0), CInteger(1), CInteger(11), CInteger(7)};
	CInteger aiv0_1_0_0_3_5[] = {CInteger(0), CInteger(1), CInteger(0), CInteger(0), CInteger(3), CInteger(5)};
	CInteger aiv0_0_1_0_9_13[] = {CInteger(0), CInteger(0), CInteger(1), CInteger(0), CInteger(9), CInteger(13)};
	CVector<CInteger> iv1_0_0_1_11_7, iv0_1_0_0_3_5, iv0_0_1_0_9_13;

	CInteger aiv1_0_0__1_78_9[] = {CInteger(1), CInteger(0), CInteger(0), CInteger(-1), CInteger(78), CInteger(9)};
	CInteger aiv0_1__1_0_17__6[] = {CInteger(0), CInteger(1), CInteger(-1), CInteger(0), CInteger(17), CInteger(-6)};
	CVector<CInteger> iv1_0_0__1_78_9, iv0_1__1_0_17__6;

	CInteger aiv0_0_0__1_87_22[] = {CInteger(0), CInteger(0), CInteger(0), CInteger(-1), CInteger(87), CInteger(22)};
	CInteger aiv1_0__1_0_28_1[] = {CInteger(1), CInteger(0), CInteger(-1), CInteger(0), CInteger(28), CInteger(1)};
	CVector<CInteger> iv0_0_0__1_87_22, iv1_0__1_0_28_1;

	CInteger aiv0_1_1_1_7[] = {CInteger(0), CInteger(1), CInteger(1), CInteger(1), CInteger(7)};

	CInteger aiv1_0_0__1_3[] = {CInteger(1), CInteger(0), CInteger(0), CInteger(-1), CInteger(3)};
	CInteger aiv0_1__1_0_2[] = {CInteger(0), CInteger(1), CInteger(-1), CInteger(0), CInteger(2)};

	CInteger aiv0_1__1__1_12[] = {CInteger(0), CInteger(1), CInteger(-1), CInteger(-1), CInteger(12)};
	CVector<CInteger> iv0_1_1_1_7, iv1_0_0__1_3, iv0_1__1_0_2, iv0_1__1__1_12;

	void test_ModuleRelation_initializeIntegerVectors() {
		iv1_0_0_1 = CVector<CInteger>(4, aiv1_0_0_1);
		iv0_1_0_0 = CVector<CInteger>(4, aiv0_1_0_0);
		iv0_0_1_0 = CVector<CInteger>(4, aiv0_0_1_0);

		iv1_0_0_0_0_1 = CVector<CInteger>(6, aiv1_0_0_0_0_1);
		iv0_1_0_0_0_0 = CVector<CInteger>(6, aiv0_1_0_0_0_0);
		iv0_0_1_0_0_0 = CVector<CInteger>(6, aiv0_0_1_0_0_0);

		iv1_0_0__1 = CVector<CInteger>(4, aiv1_0_0__1);
		iv0_1__1_0 = CVector<CInteger>(4, aiv0_1__1_0);

		iv0_0_0__1 = CVector<CInteger>(4, aiv0_0_0__1);
		iv1_0__1_0 = CVector<CInteger>(4, aiv1_0__1_0);

		iv0_1_1_1 = CVector<CInteger>(4, aiv0_1_1_1);
		iv0_1__1__1 = CVector<CInteger>(4, aiv0_1__1__1);

		iv1_0_0_1_11_7 = CVector<CInteger>(6, aiv1_0_0_1_11_7);
		iv0_1_0_0_3_5 = CVector<CInteger>(6, aiv0_1_0_0_3_5);
		iv0_0_1_0_9_13 = CVector<CInteger>(6, aiv0_0_1_0_9_13);

		iv1_0_0__1_78_9 = CVector<CInteger>(6, aiv1_0_0__1_78_9);
		iv0_1__1_0_17__6 = CVector<CInteger>(6, aiv0_1__1_0_17__6);

		iv0_0_0__1_87_22 = CVector<CInteger>(6, aiv0_0_0__1_87_22);
		iv1_0__1_0_28_1 = CVector<CInteger>(6, aiv1_0__1_0_28_1);

		iv0_1_1_1_7 = CVector<CInteger>(5, aiv0_1_1_1_7);

		iv1_0_0__1_3 = CVector<CInteger>(5, aiv1_0_0__1_3);
		iv0_1__1_0_2 = CVector<CInteger>(5, aiv0_1__1_0_2);

		iv0_1__1__1_12 = CVector<CInteger>(5, aiv0_1__1__1_12);
	}


	void CModuleRelation_TransitionOperation_1() {
		IModuleRelation a = IMR_Con(3, 1, 0, 3, IV_(4,  1,  0,  0, -1),
												IV_(4,  0,  1,  0,  0),
												IV_(4,  0,  0,  1,  0));

		IModuleRelation b = IMR_Con(1, 3, 0, 1, IV_(4,  1,  0,  0,  1));

		IModuleRelation exp = IMR_Con(3, 3, 0, 3, IV_(6,  1,  0,  0,  0,  0,  1),
				                                  IV_(6,  0,  1,  0,  0,  0,  0),
				                                  IV_(6,  0,  0,  1,  0,  0,  0));

		transitionOperation(a, b).print();
		exp.print();
		ASSERT(exp == transitionOperation(a, b));
	}

	void CModuleRelation_TransitionOperation_2() {
		test_ModuleRelation_initializeIntegerVectors();

		IModuleRelation a = IMR_Con(2, 2, 0, 3
		, IV_(4,  1,  0,  0, -1)
		, IV_(iv0_1_0_0)
		, IV_(4,  0,  0, -1,  0));

		IModuleRelation b = IMR_Con(2, 2, 0, 2
		, IV_(iv1_0_0__1)
		, IV_(iv0_1__1_0));

		IModuleRelation exp = IMR_Con(2, 2, 0, 3
		, IV_(iv0_0_0__1)
		, IV_(iv1_0__1_0)
		, IV_(iv0_1_0_0));

		ASSERT(exp == transitionOperation(a, b));
	}

	void CModuleRelation_TransitionOperation_3() {
		test_ModuleRelation_initializeIntegerVectors();

		IModuleRelation a = IMR_Con(2, 2, 0, 3
		, IV_(4,  1,  0,  0, -1)
		, IV_(iv0_1_0_0)
		, IV_(4,  0,  0, -1,  0));

		IModuleRelation b = IMR_Con(2, 2, 0, 2
		, IV_(iv1_0_0__1)
		, IV_(iv0_1__1_0));

		IModuleRelation exp = IMR_Con(2, 2, 0, 3
		, IV_(iv0_0_0__1)
		, IV_(iv1_0__1_0)
		, IV_(iv0_1_0_0));

		ASSERT(exp == transitionOperation(a, b));
	}

	void CModuleRelation_TransitionOperation_4() {
		test_ModuleRelation_initializeIntegerVectors();

		IModuleRelation a = IMR_Con(2, 2, 0, 1, IV_(4,  0,  1, -1, -1));

		IModuleRelation b = IMR_Con(2, 2, 0, 2, IV_(iv1_0_0__1)
		, IV_(iv0_1__1_0));

		IModuleRelation exp = IMR_Con(2, 2, 0, 1, IV_(iv0_1__1__1));

		ASSERT(exp == transitionOperation(a, b));
	}

	void CModuleRelation_TransitionOperation_5() {
		test_ModuleRelation_initializeIntegerVectors();

		IModuleRelation a = IMR_Con(2, 2, 2, 3, IV_(6, 11,  7,  1,  0,  0, -1)
		, IV_(6,  3,  5,  0,  1,  0,  0)
		, IV_(6,  9, 13,  0,  0, -1,  0));

		IModuleRelation b = IMR_Con(2, 2, 2, 2, IV_(6, 78,  9,  1,  0,  0, -1),
				                                IV_(6, 17, -6,  0,  1, -1,  0));

		IModuleRelation exp = IMR_Con(2, 2, 2, 3, IV_(6, 87, 22,  0,  0,  0, -1),
				                                  IV_(6, 28,  1,  1,  0, -1,  0),
				                                  IV_(6,  3,  5,  0,  1,  0,  0));

		ASSERT(exp == transitionOperation(a, b));
	}

	void CModuleRelation_TransitionOperation_6() {
		test_ModuleRelation_initializeIntegerVectors();

		IModuleRelation a = IMR_Con(2, 2, 1, 1, IV_(5,  7,  0,  1, -1, -1));

		IModuleRelation b = IMR_Con(2, 2, 1, 2, IV_(5,  3,  1,  0,  0, -1),
				                                IV_(5,  2,  0,  1, -1,  0));

		IModuleRelation exp = IMR_Con(2, 2, 1, 1, IV_(5, 12,  0,  1, -1, -1));

		ASSERT(exp == transitionOperation(a, b));
	}

	void CModuleRelation_ReflexiveTransitiveClosure_A() {
		IModuleRelation a = IMR_Con(3, 3, 1, 2, IV_(7,  0,  1,  0,  0,  0, -1,  0),
											    IV_(7,  0,  0,  1,  0,  0,  0, -1));
		IModule exp = IMod_Con(0, 3);

		ASSERT_EQ(exp, a.getReflexiveTransitiveClosure<IModule>());
	}

	void CModuleRelation_ReflexiveTransitiveClosure_B() {
		IModuleRelation a = IMR_Con(3, 3, 1, 3, IV_(7,  0,  1,  0,  0,  0, -1,  0),
											    IV_(7,  0,  0,  1,  0, -1,  0,  0),
											    IV_(7,  0,  0,  0,  1,  0,  0, -1));
		IModule exp = IMod_Con(2, IV_(3,  0,  0,  1),
				                  IV_(3,  1,  1,  0));

		ASSERT_EQ(exp, a.getReflexiveTransitiveClosure<IModule>());
	}

	void CModuleRelation_ReflexiveTransitiveClosure_C() {
		IModuleRelation a = IMR_Con(3, 3, 1, 3, IV_(7,  0,  1,  0,  0,  0, -1,  0),
											    IV_(7,  0,  0,  1,  0,  0,  0, -1),
											    IV_(7,  0,  0,  0,  1, -1,  0,  0));
		IModule exp = IMod_Con(1, IV_(3,  1,  1,  1));

		ASSERT(exp == a.getReflexiveTransitiveClosure<IModule>());
	}

	void CModuleRelation_ReflexiveTransitiveClosure_D() {
		int d = 25;
		IModuleRelation a = IModuleRelation(d, d, CFlat<IModule>(1, IModule(C, matrixToDescriptorSet(
				IMatrix::getUnitMatrix(d)
	    		<< -IMatrix::getUnitMatrix(d).getSubMatrixColumns(d - 1, 1)
			    << -IMatrix::getUnitMatrix(d).getSubMatrixColumns(0, d - 1)
			    << IMatrix::getZeroMatrix(d, 1)
	    ))));
		//IModuleRelation exp = IMR_Con(3, 3, 1, 1, IV_(7,  1,  1,  1, -1, -1, -1, 0));

//		a.calculateReflexiveTransitiveClosure().print();

		//(!a.calculateReflexiveTransitiveClosure()).print();

//		ASSERT(exp == a.calculateReflexiveTransitiveClosure());
	}
}

cute::suite* Test_ModuleRelation_runSuite(){
	cute::suite &s = *(new cute::suite("ModuleRelation"));

	s.push_back(CUTE(Test_ModuleRelation::CModuleRelation_TransitionOperation_1));
	s.push_back(CUTE(Test_ModuleRelation::CModuleRelation_TransitionOperation_2));
	s.push_back(CUTE(Test_ModuleRelation::CModuleRelation_TransitionOperation_3));
	s.push_back(CUTE(Test_ModuleRelation::CModuleRelation_TransitionOperation_4));
	s.push_back(CUTE(Test_ModuleRelation::CModuleRelation_TransitionOperation_5));
	s.push_back(CUTE(Test_ModuleRelation::CModuleRelation_TransitionOperation_6));

	s.push_back(CUTE(Test_ModuleRelation::CModuleRelation_ReflexiveTransitiveClosure_A));
	s.push_back(CUTE(Test_ModuleRelation::CModuleRelation_ReflexiveTransitiveClosure_B));
	s.push_back(CUTE(Test_ModuleRelation::CModuleRelation_ReflexiveTransitiveClosure_C));
	//s.push_back(CUTE(Test_ModuleRelation::CModuleRelation_ReflexiveTransitiveClosure_D));

	return &s;
}
