
#include "../cute/cute.h"
#include "../cute/ide_listener.h"

#include "Test_Funcs.h"

#include "../partyql/basicmath/scalar/Rational.h"
#include "../partyql/basicmath/scalar/Interval.h"
#include "../partyql/basicmath/Polyheder.h"

using namespace AlgoTrans;

namespace Test_FourierMotzkin {

	typedef CRational<I> IRat;
	typedef CInterval<IRat> IRatIntv;

	IRatIntv IRatIntV(int lowerBound, int upperBound) {
		return IRatIntv(IRat(lowerBound), IRat(upperBound));
	}

	IRatIntv IRatIntV(bool lowerBounded, int bound) {
		return IRatIntv(lowerBounded, IRat(bound));
	}
/*
	void A() {
		IPolyheder pIn = IP(C, false, 4, IV_(3,  2,  1, -1),  // x - y >= -2
				               IV_(3,  2, -1,  1),  // x - y <=  2
				               IV_(3,  8, -1, -1),  // x + y <=  8
				               IV_(3, -4,  1,  1)); // x + y >=  4

		IRatIntv pProjX = pIn.getProjectionFourierMotzkin(0).asInterval();
		IRatIntv pProjY = pIn.getProjectionFourierMotzkin(1).asInterval();

		ASSERT(pProjX == IRatIntV(1, 5));
		ASSERT(pProjY == IRatIntV(1, 5));

		IRatIntv pProjX_Jones = pIn.getProjectionJones(0).asInterval();
		IRatIntv pProjY_Jones = pIn.getProjectionJones(1).asInterval();

		ASSERT(pProjX == pProjX_Jones);
		ASSERT(pProjY == pProjY_Jones);
	}

	void B() { // A with additional redundant constraints
		IPolyheder pIn = IP(C, false, 6, IV_(3,  4,  2, -2),  // 2x - 2y >= -4
				 			   IV_(3,  2,  1, -1),  // x - y >= -2
				               IV_(3,  2, -1,  1),  // x - y <=  2
				               IV_(3,  8, -1, -1),  // x + y <=  8
				               IV_(3, 20, -1, -1),  // x + y <= 20
				               IV_(3, -4,  1,  1)); // x + y >=  4

		IRatIntv pProjX = pIn.getProjectionFourierMotzkin(0).asInterval();
		IRatIntv pProjY = pIn.getProjectionFourierMotzkin(1).asInterval();

		ASSERT(pProjX == IRatIntV(1, 5));
		ASSERT(pProjY == IRatIntV(1, 5));

		IRatIntv pProjX_Jones = pIn.getProjectionJones(0).asInterval();
		IRatIntv pProjY_Jones = pIn.getProjectionJones(1).asInterval();

		ASSERT(pProjX == pProjX_Jones);
		ASSERT(pProjY == pProjY_Jones);
	}

	void C() { // A with one constraint removed such that domain is unbounded in one direction
		IPolyheder pIn = IP(C, false, 3, IV_(3,  2,  1, -1),  // x - y >= -2
				               IV_(3,  2, -1,  1),  // x - y <=  2
				               IV_(3, -4,  1,  1)); // x + y >=  4

		IRatIntv pProjX = pIn.getProjectionFourierMotzkin(0).asInterval();
		IRatIntv pProjY = pIn.getProjectionFourierMotzkin(1).asInterval();

		ASSERT(pProjX == IRatIntV(true, 1)); // [1, +inf[
		ASSERT(pProjY == IRatIntV(true, 1)); // [1, +inf[

		IRatIntv pProjX_Jones = pIn.getProjectionJones(0).asInterval();
		IRatIntv pProjY_Jones = pIn.getProjectionJones(1).asInterval();

		ASSERT(pProjX == pProjX_Jones);
		ASSERT(pProjY == pProjY_Jones);
	}*/
}


cute::suite* Test_FourierMotzkin_runSuite(){
	cute::suite& s = *(new cute::suite("Polyhedral Projection"));

	//s.push_back(CUTE(Test_FourierMotzkin::A));
	//s.push_back(CUTE(Test_FourierMotzkin::B));
	//s.push_back(CUTE(Test_FourierMotzkin::C));

	return &s;
}


