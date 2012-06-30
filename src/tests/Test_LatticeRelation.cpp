#include "../cute/cute.h"
#include "../cute/ide_listener.h"
#include "../cute/cute_suite.h"

#include "Test_Funcs.h"

#include "../partyql/graph/Graph.h"
#include "../partyql/basicmath/scalar/Integer.h"
#include "../partyql/basicmath/Vektor.h"
#include "../partyql/basicmath/Lattice.h"
#include "../partyql/basicmath/SetRelation.h"

using namespace AlgoTrans;

namespace Test_LatticeRelation {
	void CLatticeRelation_TransitionOperation_1() {
		ILatticeRelation a = ILR(1, 1, 1, 2, IV_(3,  1,  1,  0),
										     IV_(3,  0,  0,  1));

		ILatticeRelation b = ILR(1, 1, 1, 1, IV_(3,  1,  1,  1));

		ILatticeRelation exp = ILR(1, 1, 1, 1, IV_(3,  1,  1,  1));

		ASSERT(exp == transitionOperation(a, b));
	}

	void CLatticeRelation_TransitionOperation_2() {
		ILatticeRelation a = ILR(1, 1, 2, 3, IV_(4,  1,  1,  0,  0),
										     IV_(4,  0,  0,  1,  0),
										     IV_(4,  0,  0,  0,  1));

		ILatticeRelation b = ILR(1, 1, 2, 2, IV_(4,  1,  1,  0,  1),
										     IV_(4,  0,  0,  1,  0));

		ILatticeRelation exp = ILR(1, 1, 2, 2, IV_(4,  1,  1,  0,  1),
										       IV_(4,  0,  0,  1,  0));

		ASSERT(exp == transitionOperation(a, b));
	}

	void CLatticeRelation_TransitionOperation_3() {
		ILatticeRelation a = ILR(2, 2, 2, 4, IV_(6,  1,  0,  1,  0,  0,  0),
										     IV_(6,  0,  1,  0,  1,  0,  0),
										     IV_(6,  0,  0,  0,  0,  1,  0),
										     IV_(6,  0,  0,  0,  0,  0,  1));

		ILatticeRelation b = ILR(2, 2, 2, 3, IV_(6,  1,  0,  1,  0,  0,  0),
										     IV_(6,  0,  1,  0,  1,  0,  1),
										     IV_(6,  0,  0,  0,  0,  1,  0));

		ILatticeRelation exp = ILR(2, 2, 2, 3, IV_(6,  1,  0,  1,  0,  0,  0),
							  			       IV_(6,  0,  1,  0,  1,  0,  1),
										       IV_(6,  0,  0,  0,  0,  1,  0));

		ASSERT(exp == transitionOperation(a, b));
	}
}

cute::suite* Test_LatticeRelation_runSuite(){
	cute::suite &s = *(new cute::suite("Lattice Relation"));

	s.push_back(CUTE(Test_LatticeRelation::CLatticeRelation_TransitionOperation_1));
	s.push_back(CUTE(Test_LatticeRelation::CLatticeRelation_TransitionOperation_2));
	s.push_back(CUTE(Test_LatticeRelation::CLatticeRelation_TransitionOperation_3));

	return &s;
}
