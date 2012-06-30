
#include "../cute/cute.h"
#include "../cute/ide_listener.h"

#include "Test_Funcs.h"

#include "../partyql/basicmath/scalar/Rational.h"
#include "../partyql/basicmath/scalar/Interval.h"
#include "../partyql/basicmath/Polyheder.h"
#include "../partyql/basicmath/IntegerProgramming.h"

using namespace AlgoTrans;

namespace Test_IntegerProgramming {
	void C() {
		CIpProblem<I> ip = CIpProblem<I>(2);

		ip.setDomainPolyheder(CFlat<ICone>(1, ICone_Con(4, IV_(3,  2,  1, -1),  // x - y >= -2
                                             IV_(3,  2, -1,  1),  // x - y <=  2
	                                         IV_(3,  8, -1, -1),  // x + y <=  8
	                                         IV_(3, -4,  1,  1)))); // x + y >=  4

		CVector<I::RationalType>* sol = ip.solve();

		sol->print();

		//ASSERT(affineHull == exp);
	}

	void D() {
		CPipProblem<I> ip = CPipProblem<I>(0, 2);

		ip.domConstrPolyheder = CFlat<ICone>(1, ICone_Con(4, IV_(3,  2,  1, -1),  // x - y >= -2
                                             IV_(3,  2, -1,  1),  // x - y <=  2
	                                         IV_(3,  8, -1, -1),  // x + y <=  8
	                                         IV_(3, -4,  1,  1),
	                                         IV_(3,  1,  0,  0))); // x + y >=  4

		CPipQuast<I>* sol = ip.solve();

		sol->print();

		//ASSERT(affineHull == exp);
	}
}

cute::suite* Test_IntegerProgramming_runSuite(){
	cute::suite& s = *(new cute::suite("Integer programming"));

	s.push_back(CUTE(Test_IntegerProgramming::C));

	return &s;
}



