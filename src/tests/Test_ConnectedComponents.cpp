
#include "../cute/cute.h"
#include "../cute/ide_listener.h"

#include "Test_Funcs.h"

#include "../partyql/graph/Graph.h"
#include "../partyql/basicmath/scalar/Integer.h"

using namespace AlgoTrans;

namespace Test_ConnectedComponentFinder {
	template <class E>
	bool equalSet(vector<E>& a, vector<E>& b) {
		if (a.size() != b.size()) return false;

		vector<bool> used = vector<bool>(a.size(), false);
		for (unsigned int q = 0; q < a.size(); q++) {
			bool found = false;
			for (unsigned int r = 0; r < b.size(); r++) if (!used[r]) {
				if (a[q] == b[r]) {
					used[r] = true;
					found = true;
					break;
				}
			}
			if (!found) return false;
		}

		return true;
	}

	template <class E>
	bool equalSetOfSets(vector<vector<E> >& a, vector<vector<E> >& b) {
		if (a.size() != b.size()) return false;

		vector<bool> used = vector<bool>(a.size(), false);
		for (unsigned int q = 0; q < a.size(); q++) {
			bool found = false;
			for (unsigned int r = 0; r < b.size(); r++) if (!used[r]) {
				if (equalSet<E>(a[q], b[r])) {
					used[r] = true;
					found = true;
					break;
				}
			}
			if (!found) return false;
		}

		return true;
	}

	void printComponents(vector<vector<CVertex<int, bool>* > > c) {
		printf("[");
		for (unsigned int q = 0; q < c.size(); q++) {
			if (q > 0) printf(",\n");
			printf("[");
			for (unsigned int r = 0; r < c[q].size(); r++) {
				if (r > 0) printf(", ");
				printf(CInteger(c[q][r]->getData()).toString().c_str());
			}
			printf("]");
		}
		printf("]\n");
	}

	bool equalSetOfVertexSets(vector<vector<CVertex<int, bool>* > >& a, vector<vector<CVertex<int, bool>* > >& b) {
		return equalSetOfSets<CVertex<int, bool>* >(a, b);
	}

	void C() {
		CGraph<int, bool> g;

		int vertCount = 2;
		for (int q = 0; q < vertCount; q++) g.addVertex(q);

		g.addEdge(g(0), g(0), true);

		vector<vector<CVertex<int, bool>* > > components = g.getConnectedComponents();

		// Exp
		vector<vector<CVertex<int, bool>* > > exp = vector<vector<CVertex<int, bool>* > >();

		vector<CVertex<int, bool>* > expA = vector<CVertex<int, bool>* >();
		expA.push_back(&g(0));

		vector<CVertex<int, bool>* > expB = vector<CVertex<int, bool>* >();
		expB.push_back(&g(1));

		exp.push_back(expA);
		exp.push_back(expB);

		ASSERT(equalSetOfVertexSets(exp, components));
	}

	void B() {
		CGraph<int, bool> g;

		int vertCount = 6;
		for (int q = 0; q < vertCount; q++) g.addVertex(q);

		g.addEdge(g(0), g(1), true);
		g.addEdge(g(0), g(2), true);

		g.addEdge(g(3), g(4), true);

		g.addEdge(g(5), g(5), true);

		vector<vector<CVertex<int, bool>* > > components = g.getConnectedComponents();

		// Exp
		vector<vector<CVertex<int, bool>* > > exp = vector<vector<CVertex<int, bool>* > >();

		vector<CVertex<int, bool>* > expA = vector<CVertex<int, bool>* >();
		expA.push_back(&g(0));
		expA.push_back(&g(1));
		expA.push_back(&g(2));

		vector<CVertex<int, bool>* > expB = vector<CVertex<int, bool>* >();
		expB.push_back(&g(3));
		expB.push_back(&g(4));

		vector<CVertex<int, bool>* > expC = vector<CVertex<int, bool>* >();
		expC.push_back(&g(5));

		exp.push_back(expA);
		exp.push_back(expB);
		exp.push_back(expC);

		ASSERT(equalSetOfVertexSets(exp, components));
	}

	void A() {
		CGraph<int, bool> g;

		int vertCount = 7;
		for (int q = 0; q < vertCount; q++) g.addVertex(q);

		g.addEdge(g(0), g(1), true);
		g.addEdge(g(0), g(2), true);

		g.addEdge(g(3), g(4), true);

		g.addEdge(g(5), g(5), true);

		vector<vector<CVertex<int, bool>* > > components = g.getConnectedComponents();

		// Exp
		vector<vector<CVertex<int, bool>* > > exp = vector<vector<CVertex<int, bool>* > >();

		vector<CVertex<int, bool>* > expA = vector<CVertex<int, bool>* >();
		expA.push_back(&g(0));
		expA.push_back(&g(1));
		expA.push_back(&g(2));

		vector<CVertex<int, bool>* > expB = vector<CVertex<int, bool>* >();
		expB.push_back(&g(3));
		expB.push_back(&g(4));

		vector<CVertex<int, bool>* > expC = vector<CVertex<int, bool>* >();
		expC.push_back(&g(5));

		vector<CVertex<int, bool>* > expD = vector<CVertex<int, bool>* >();
		expD.push_back(&g(6));

		exp.push_back(expA);
		exp.push_back(expB);
		exp.push_back(expC);
		exp.push_back(expD);

		ASSERT(equalSetOfVertexSets(exp, components));
	}
}

cute::suite* Test_ConnectedComponentFinder_runSuite(){
	cute::suite& s = *(new cute::suite("Connected Components Finder"));

	s.push_back(CUTE(Test_ConnectedComponentFinder::C));
	s.push_back(CUTE(Test_ConnectedComponentFinder::B));
	s.push_back(CUTE(Test_ConnectedComponentFinder::A));

	return &s;
}


