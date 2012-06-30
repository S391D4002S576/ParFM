
#include "../cute/cute.h"
#include "../cute/ide_listener.h"

#include "Test_Funcs.h"

#include "../partyql/basicmath/scalar/Rational.h"
#include "../partyql/basicmath/scalar/Interval.h"
#include "../partyql/basicmath/Polyheder.h"

using namespace AlgoTrans;

namespace Test_Cone {
	void D() { // dims: a, b, c, d
		ICone a = ICone_Con(2, IV_(5, -1, -1,  0,  1,  0),
				               IV_(5,  1,  0,  0,  0,  0)); // c >= a + 1 == c > a

		ICone b = ICone_Con(4, IV_(5,  0, -1,  0,  1,  0), // c >= a
				 			   IV_(5,  0,  1,  0, -1,  0), // c <= a
				 			   IV_(5, -1,  0, -1,  0,  1),
				               IV_(5,  1,  0,  0,  0,  0)); // d >= b + 1 == d > b

		ICone aD = !a;
		ICone bD = !b;

		//aD.print();
		//bD.print();

		//(a + b).print();
		//ASSERT((a && b) == a);

		//ASSERT((a + b) == b);
	}

	void A() {
		ICone a = ICone_Con(2, IV_(2,  2,  1),
				               IV_(2,  1,  0));  // x >= -2 , 1 >= 0

		ICone b = ICone_Con(2, IV_(2,  3,  1),
	                           IV_(2,  1,  0));  // x >= -3

		//ASSERT((a && b) == a);

		//ASSERT((a + b) == b);
	}

	void B() {
		ICone a = ICone_Con(2, IV_(2,  2,  1),
							   IV_(2, -3, -1),
							   IV_(2,  1,  0));

		//(!a).print();

		//ASSERT((a && b) == a);

		//ASSERT((a + b) == b);
	}

	void C() {
		ICone a = ICone_Con(5, IV_(3,  2,  1, -1),  // x - y >= -2
			            IV_(3,  2, -1,  1),  // x - y <=  2
				        IV_(3,  8, -1, -1),  // x + y <=  8
				        IV_(3, -4,  1,  1),
				        IV_(3,  1,  0,  0)); // x + y >=  4

		ICone b = ICone_Con(5, IV_(3,  0,  1, -1),  // x - y >=  0
			            IV_(3,  4, -1,  1),  // x - y <=  4
				        IV_(3, 10, -1, -1),  // x + y <= 10
				        IV_(3, -6,  1,  1),
				        IV_(3,  1,  0,  0)); // x + y >=  6

		//(a && b).print();

		//(a + b).print();
	}
}


cute::suite* Test_Cone_runSuite(){
	cute::suite& s = *(new cute::suite("Cone"));

	s.push_back(CUTE(Test_Cone::D));
	s.push_back(CUTE(Test_Cone::A));
	s.push_back(CUTE(Test_Cone::B));
	s.push_back(CUTE(Test_Cone::C));

	return &s;
}


