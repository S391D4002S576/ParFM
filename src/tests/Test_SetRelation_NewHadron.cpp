#include "../cute/cute.h"
#include "../cute/ide_listener.h"
#include "../cute/cute_suite.h"

#include "Test_Funcs.h"

#include "../partyql/graph/Graph.h"
#include "../partyql/core/scalar/Integer.h"
#include "../partyql/core/SetRelation.h"

using namespace AlgoTrans;

namespace Test_SetRelation {
	// From Lattice Relation
	void SetRelation_TransitionOperation_1() { // Lattice relation test 1
		ISetRelation a = ISR(true, G, 0, 1, 1, 2, IV_(3,  1,  1,  0),
  									     IV_(3,  0,  0,  1));
		ISetRelation b = ISR(true, G, 0, 1, 1, 1, IV_(3,  1,  1,  1));
		ISetRelation exp = ISR(true, G, 0, 1, 1, 1, IV_(3,  1,  1,  1));

		ISetRelation::Conjunctive_CombinationOperation combOp = ISetRelation::Conjunctive_CombinationOperation();
		ISetRelation act = transitionOperation(combOp, a, b);
		//exp.print();
		//act.print();

		ASSERT_EQ(exp, act);
	}

	void SetRelation_TransitionOperation_2() { // Lattice relation test 2
		ISetRelation a = ISR(true, G, 1, 1, 1, 3, IV_(4,  1,  1,  0,  0),
									        IV_(4,  0,  0,  1,  0),
									        IV_(4,  0,  0,  0,  1));

		ISetRelation b = ISR(true, G, 1, 1, 1, 2, IV_(4,  1,  1,  0,  1),
									        IV_(4,  0,  0,  1,  0));

		ISetRelation exp = ISR(true, G, 1, 1, 1, 2, IV_(4,  1,  1,  0,  1),
									          IV_(4,  0,  0,  1,  0));

		ISetRelation::Conjunctive_CombinationOperation combOp = ISetRelation::Conjunctive_CombinationOperation();
		ISetRelation act = transitionOperation(combOp, a, b);
		//act.print();

		ASSERT_EQ(exp, act);
	}

	void SetRelation_TransitionOperation_3() { // Lattice relation test 3
		ISetRelation a = ISR(true, G, 1, 2, 2, 4, IV_(6,  1,  0,  1,  0,  0,  0),
										     IV_(6,  0,  1,  0,  1,  0,  0),
										     IV_(6,  0,  0,  0,  0,  1,  0),
										     IV_(6,  0,  0,  0,  0,  0,  1));

		ISetRelation b = ISR(true, G, 1, 2, 2, 3, IV_(6,  1,  0,  1,  0,  0,  0),
										     IV_(6,  0,  1,  0,  1,  0,  1),
										     IV_(6,  0,  0,  0,  0,  1,  0));

		ISetRelation exp = ISR(true, G, 1, 2, 2, 3, IV_(6,  1,  0,  1,  0,  0,  0),
												IV_(6,  0,  1,  0,  1,  0,  1),
										       IV_(6,  0,  0,  0,  0,  1,  0));

		ISetRelation::Conjunctive_CombinationOperation combOp = ISetRelation::Conjunctive_CombinationOperation();
		ISetRelation act = transitionOperation(combOp, a, b);
		//act.print();

		ASSERT_EQ(exp, act);
	}

	void SetRelation_TransitionOperation_4() { // Module relation test 1
		ISetRelation a = ISR(true, C, 0, 3, 1, 3, IV_(5,  0,  1,  0,  0, -1),
												IV_(5,  0,  0,  1,  0,  0),
												IV_(5,  0,  0,  0,  1,  0));

		ISetRelation b = ISR(true, C, 0, 1, 3, 1, IV_(5,  0,  1,  0,  0,  1));

		ISetRelation exp = ISR(true, C, 0, 3, 3, 3, IV_(7,  0,  1,  0,  0,  0,  0,  1),
				                                  IV_(7,  0,  0,  1,  0,  0,  0,  0),
				                                  IV_(7,  0,  0,  0,  1,  0,  0,  0));

		ISetRelation::Conjunctive_CombinationOperation combOp = ISetRelation::Conjunctive_CombinationOperation();
		ISetRelation act = transitionOperation(combOp, a, b);
		//act.print();

		ASSERT_EQ(exp, act);
	}

	void SetRelation_TransitionOperation_5() { // Module relation test 1
		ISetRelation a = ISR(true, C, 0, 1, 1, 1,  IV_(3,  1,  1,  -1)) *
				         ISR(false, C, 0, 1, 1, 2, IV_(3,  0,  1,  0),
						 						   IV_(3,  1,  0,  0));

		ISetRelation exp = ISR(true, C, 0, 1, 1, 1,  IV_(3,  2,  1,  -1)) *
				           ISR(false, C, 0, 1, 1, 2, IV_(3,  0,  1,  0),
						 						     IV_(3,  1,  0,  0));

		ISetRelation::Conjunctive_CombinationOperation combOp = ISetRelation::Conjunctive_CombinationOperation();
		ISetRelation act = transitionOperation(combOp, a, a);
		//act.print();
		//exp.print();

		ASSERT_EQ(exp, act);
	}

	/*
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
	 */

	void SetRelation_TransitiveClosure_A() {
		ISetRelation a = ISR(true, C, 0, 3, 3, 2, IV_(7,  0,  1,  0,  0,  0, -1,  0),
										          IV_(7,  0,  0,  1,  0,  0,  0, -1));
		ISetRelation exp = ISR(true, C, 0, 3, 3, 0);

		ISetRelation act = a.getTransitiveClosure();
		//act.print();
		//act.print();

		ASSERT_EQ(exp, act);
	}

	void SetRelation_TransitiveClosure_B() {
		ISetRelation a = ISR(true, C, 0, 3, 3, 3, IV_(7,  0,  1,  0,  0,  0, -1,  0),
										          IV_(7,  0,  0,  1,  0, -1,  0,  0),
										          IV_(7,  0,  0,  0,  1,  0,  0, -1));
		ISetRelation exp = ISR(true, C, 0, 3, 3, 2, IV_(7,  0,  1,  1,  0, -1, -1,  0),
											  IV_(7,  0,  0,  0,  1,  0,  0, -1));

		ISetRelation act = a.getTransitiveClosure();
		//act.print();

		ASSERT_EQ(exp, act);
	}

	void SetRelation_TransitiveClosure_C() {
		ISetRelation a = ISR(true, C, 0, 3, 3, 3, IV_(7,  0,  1,  0,  0,  0, -1,  0),
										    IV_(7,  0,  0,  1,  0,  0,  0, -1),
										    IV_(7,  0,  0,  0,  1, -1,  0,  0));
		ISetRelation exp = ISR(true, C, 0, 3, 3, 1, IV_(7,  0,  1,  1,  1,  -1,  -1,  -1));

		ISetRelation act = a.getTransitiveClosure();
		//act.print();

		ASSERT_EQ(exp, act);
	}

	void SetRelation_TransitiveClosure_D() {
		ISetRelation a = ISR(false, C, 1, 1, 1, 6, IV_(4, -1,  0, -1,  1),
		                                           IV_(4,  1,  0,  0,  0),
										           IV_(4,  0,  0,  0,  1),
										           IV_(4,  0,  0,  1,  0),
										           IV_(4,  0,  1, -1,  0),
										           IV_(4,  0,  1,  0, -1));
		ISetRelation exp = ISR(false, C, 1, 1, 1, 1, IV_(4,  0,  0,  -1,  1));

		ISetRelation act = a.getTransitiveClosure();
		//act.print(true);

		ASSERT_EQ(exp, act);
	}

	void CModuleRelation_ReflexiveTransitiveClosure_D() {
		int d = 2;
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

cute::suite* Test_SetRelation_NewHadron_runSuite(){
	cute::suite &s = *(new cute::suite("Set Relation"));

	s.push_back(CUTE(Test_SetRelation::SetRelation_TransitionOperation_4));

	s.push_back(CUTE(Test_SetRelation::SetRelation_TransitionOperation_1));
	s.push_back(CUTE(Test_SetRelation::SetRelation_TransitionOperation_2));
	s.push_back(CUTE(Test_SetRelation::SetRelation_TransitionOperation_3));

	s.push_back(CUTE(Test_SetRelation::SetRelation_TransitionOperation_5));

	s.push_back(CUTE(Test_SetRelation::SetRelation_TransitiveClosure_A));
	s.push_back(CUTE(Test_SetRelation::SetRelation_TransitiveClosure_B));
	s.push_back(CUTE(Test_SetRelation::SetRelation_TransitiveClosure_C));

	s.push_back(CUTE(Test_SetRelation::SetRelation_TransitiveClosure_D));

	return &s;
}
