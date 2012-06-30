
#include "../cute/cute.h"
#include "../cute/ide_listener.h"

#include "Test_Funcs.h"

#include "../partyql/graph/Graph.h"

using namespace AlgoTrans;

namespace Test_StronglyConnectedComponentFinder {
	void A() {
		//CGraph<>
	}
}

cute::suite* Test_StronlgyConnectedComponentFinder_runSuite(){
	cute::suite& s = *(new cute::suite("Strongly Connected Components Finder"));

	s.push_back(CUTE(Test_StronglyConnectedComponentFinder::A));

	return &s;
}


