/*
#include "cute.h"
#include "ide_listener.h"

#include "Test_Funcs.h"

#include "../partyql/basicmath/FinitePartition.h"

using namespace AlgoTrans;

namespace Test_FinitePartittion {
	void coarsen_small() {
		CFiniteSet a = CFiniteSet(7);

		vector<int> xMap;
		xMap.push_back(0);
		xMap.push_back(1);
		xMap.push_back(2);
		xMap.push_back(2);
		xMap.push_back(3);
		xMap.push_back(4);
		xMap.push_back(5);
		CFinitePartition x = CFinitePartition(a, xMap);

		vector<int> yMap;
		yMap.push_back(3);
		yMap.push_back(2);
		yMap.push_back(0);
		yMap.push_back(1);
		yMap.push_back(3);
		yMap.push_back(4);
		yMap.push_back(5);
		CFinitePartition y = CFinitePartition(a, yMap);

		vector<int> zMap;
		zMap.push_back(2);
		zMap.push_back(1);
		zMap.push_back(0);
		zMap.push_back(0);
		zMap.push_back(2);
		zMap.push_back(3);
		zMap.push_back(4);
		CFinitePartition z = CFinitePartition(a, zMap);

		vector<int> gMap;
		gMap.push_back(0);
		gMap.push_back(1);
		gMap.push_back(2);
		gMap.push_back(3);
		gMap.push_back(2);
		gMap.push_back(4);
		gMap.push_back(5);
		CFinitePartition g = CFinitePartition(a, gMap);

		vector<int> hMap;
		hMap.push_back(0);
		hMap.push_back(1);
		hMap.push_back(0);
		hMap.push_back(0);
		hMap.push_back(0);
		hMap.push_back(2);
		hMap.push_back(3);
		CFinitePartition h = CFinitePartition(a, hMap);

		// Coarsen, ||
		ASSERT(z == (x || y));
		ASSERT(z == (y || x));

		ASSERT(h == (z || g));
		ASSERT(h == (g || z));

		// Refine, &&
		ASSERT(x == ((x || y) && x));
		ASSERT(y == ((y || x) && y));
		ASSERT(x == (x && (x || y)));
		ASSERT(y == (y && (y || x)));

		ASSERT(z == ((z || g) && z));
		ASSERT(g == ((g || z) && g));

		ASSERT(x == ((x && y) || x));
		ASSERT(y == ((y && x) || y));

		ASSERT(z == ((z && g) || z));
		ASSERT(g == ((g && z) || g));

		// IsNoCoarserThan, <=
		ASSERT(x <= (x || y));
		ASSERT(y <= (x || y));
		ASSERT((x && y) <= (x || y));

		ASSERT(z <= (z || g));
		ASSERT(g <= (z || g));
		ASSERT((z && g) <= (z || g));
		ASSERT((z && g) <= z);
		ASSERT((z && g) <= g);
	}

	void coarsen_larger() {
		CFiniteSet a = CFiniteSet(3*5*7);

		vector<int> xMap, yMap, zMap;
		for (int q = 0; q < a.getElementCount(); q++) {
			xMap.push_back(q % (3*7));
			yMap.push_back(q % (7*5));
			zMap.push_back(q % 7);
		}
		CFinitePartition x = CFinitePartition(a, xMap);
		CFinitePartition y = CFinitePartition(a, yMap);
		CFinitePartition z = CFinitePartition(a, zMap);

		ASSERT(z == (x || y));
		ASSERT(z == (y || x));

		ASSERT(x == ((x || y) && x));
		ASSERT(z == ((x || y) && z));
		ASSERT(y == ((x || y) && y));

		ASSERT(x == ((x && y) || x));
		ASSERT(z == ((x && y) || z));
		ASSERT(y == ((x && y) || y));
	}
}

cute::suite* Test_FinitePartition_runSuite(){
	cute::suite& s = *(new cute::suite("Finite Partition"));

	s.push_back(CUTE(Test_FinitePartittion::coarsen_small));
	s.push_back(CUTE(Test_FinitePartittion::coarsen_larger));

	return &s;
}


*/
