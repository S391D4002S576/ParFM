#include "../cute/cute.h"
#include "../cute/ide_listener.h"
#include "../cute/cute_suite.h"

#include "Test_Funcs.h"

#include "../partyql/graph/Graph.h"
#include "../partyql/basicmath/scalar/Integer.h"
#include "../partyql/basicmath/Vektor.h"
#include "../partyql/basicmath/Module.h"
#include "../partyql/basicmath/SetRelation.h"
#include "../partyql/RADGraph.h"

using namespace AlgoTrans;

namespace Test_RADGraph {
	typedef CInteger I;
	typedef CVector<I> IntegerVector;

	typedef CVertex<CPolyheder<I>, CModuleRelation<I>::T > RADVertex;
	typedef CRADGraph<I, IModuleRelation> RADGraph;

/*	void CRADGraph_TransitiveClosure() {
		RADGraph g = RADGraph();
		vector<RADVertex* > v = vector<RADVertex* >(4);
		for (unsigned int q = 0; q < v.size(); q++) v[q] = &g.addVertex();

		IModuleRelation mA = IMR_Con(2, 2, 2, 3, IV_(6, 11,  7,  1,  0,  0, -1),
				                                 IV_(6,  3,  5,  0,  1,  0,  0),
				                                 IV_(6,  9, 13,  0,  0, -1,  0));

		IModuleRelation mB = IMR_Con(2, 2, 2, 2, IV_(6, 78,  9,  1,  0,  0,  1),
				                                 IV_(6, 17, -6,  0,  1,  1,  0));

		IModuleRelation mC = IMR_Con(2, 3, 2, 1, IV_(7,  7,  2, -1, -2,  0,  1, -5));

		// Expected result of B ° A
		IModuleRelation mAB = IMR_Con(2, 2, 2, 3, IV_(6, 87, 22,  0,  0,  0,  1),
				                                  IV_(6, 28,  1,  1,  0,  1,  0),
				                                  IV_(6,  3,  5,  0,  1,  0,  0));

		// Expected result of C ° B
		IModuleRelation mBC = IMR_Con(2, 3, 2, 1, IV_(7, 180, 14,  2,  1,  0,  1, -5));

		// Expected result of C ° B ° A
		IModuleRelation mABC = IMR_Con(2, 3, 2, 2, IV_(7, 209, 47,  1,  0,  0,  1, -5),
				                                   IV_(7,   3,  5,  0,  1,  0,  0,  0));

		g.addEdge(*v[0], *v[1], mA);
		g.addEdge(*v[1], *v[2], mB);
		g.addEdge(*v[2], *v[3], mC);

		RADGraph* radg = (RADGraph*) g.calcConjunctiveTransitiveClosure();

		ASSERT(radg->getUniqueEdgeBetween(0, 1).getData() == mA);
		ASSERT(radg->getEdgesBetween(1, 0).size() == 0);
		ASSERT(radg->getUniqueEdgeBetween(1, 2).getData() == mB);
		ASSERT(radg->getEdgesBetween(2, 1).size() == 0);
		ASSERT(radg->getUniqueEdgeBetween(0, 2).getData() == mAB);
		ASSERT(radg->getUniqueEdgeBetween(1, 3).getData() == mBC);
		ASSERT(radg->getUniqueEdgeBetween(0, 3).getData() == mABC);
	}*/

	/* Depends on HNF reduction through Module-equality
	 *
	 * for x in ... do begin
	 *   B[x] := 2*C[x] * D[x];
	 * end;
	 *
	 * for (x, y, z) in ... do begin
	 *   A[x][y][z] := B[x][y] + C[y][z] + D[z][x]
	 * end;
	 *//*
	void CRADGraph_TransitiveClosure_2() {
		RADGraph g = RADGraph();
		vector<RADVertex* > v = vector<RADVertex* >(4);
		for (unsigned int q = 0; q < v.size(); q++) v[q] = &g.addVertex();

		ofstream ofs("examples/CRADGraph_TransitiveClosure_2.xhtml");
		ofs << "<html xmlns=\"http://www.w3.org/1999/xhtml\"><body>";
		g.osHtml = &ofs;

		g.addEdge(v[1], v[1], IMR(1, 1, 1, 1, IV_(3,  1, -1,  0)));
		g.addEdge(v[2], v[2], IMR(1, 1, 1, 1, IV_(3,  1, -1,  0)));
		g.addEdge(v[3], v[3], IMR(1, 1, 1, 1, IV_(3,  1, -1,  0)));

		g.addEdge(v[1], v[2], IMR(1, 1, 1, 1, IV_(3,  1, -1,  0)));
		g.addEdge(v[2], v[1], IMR(1, 1, 1, 1, IV_(3,  1, -1,  0)));
		g.addEdge(v[1], v[3], IMR(1, 1, 1, 1, IV_(3,  1, -1,  0)));
		g.addEdge(v[3], v[1], IMR(1, 1, 1, 1, IV_(3,  1, -1,  0)));
		g.addEdge(v[2], v[3], IMR(1, 1, 1, 1, IV_(3,  1, -1,  0)));
		g.addEdge(v[3], v[2], IMR(1, 1, 1, 1, IV_(3,  1, -1,  0)));

		g.addEdge(v[1], v[0], IMR(1, 3, 1, 2, IV_(5,  1, -1,  0,  0,  0),
				                              IV_(5,  1,  0,  0, -1,  0)));
		g.addEdge(v[0], v[1], IMR(3, 1, 1, 2, IV_(5,  1,  0,  0, -1,  0),
		         							  IV_(5,  0,  0,  1, -1,  0)));
		g.addEdge(v[2], v[0], IMR(1, 3, 1, 2, IV_(5,  1,  0, -1,  0,  0),
				                              IV_(5,  1,  0,  0, -1,  0)));
		g.addEdge(v[0], v[2], IMR(3, 1, 1, 2, IV_(5,  0,  1,  0, -1,  0),
		 		                              IV_(5,  0,  0,  1, -1,  0)));
		g.addEdge(v[3], v[0], IMR(1, 3, 1, 2, IV_(5,  1, -1,  0,  0,  0),
				                              IV_(5,  1,  0, -1,  0,  0)));
		g.addEdge(v[0], v[3], IMR(3, 1, 1, 2, IV_(5,  1,  0,  0, -1,  0),
				                              IV_(5,  0,  1,  0, -1,  0)));

		g.addEdge(v[0], v[0], IMR(3, 3, 1, 3, IV_(7,  1,  0,  0, -1,  0,  0,  0),
				                              IV_(7,  0,  1,  0,  0, -1,  0,  0),
								              IV_(7,  0,  0,  1,  0,  0, -1,  0)));

		RADGraph* radg = (RADGraph*) g.calcConjunctiveTransitiveClosure();

		vector<RADVertex* > w = vector<RADVertex* >();
		for (unsigned int q = 0; q < v.size(); q++) w.push_back(&radg->getVertex(q));

		vector<vector<vector<RADEdge*> > > resultEdges = extractResultEdges(radg, w);

		ofs << "</body></html>"; ofs.close();
	}*/

	/* Depends on HNF reduction through Module-equality
	 *
	 * for x in ... do begin
	 *   B[x] := 2*C[x];
	 * end;
	 *
	 * for (x, y) in ... do begin
	 *   A[x][y] := B[x] + C[y]
	 * end;
	 *//*
	void CRADGraph_TransitiveClosure_3() {
		RADGraph g = RADGraph();
		vector<RADVertex* > v = vector<RADVertex* >(3);
		for (unsigned int q = 0; q < v.size(); q++) v[q] = &g.addVertex();

		ofstream ofs("examples/CRADGraph_TransitiveClosure_3.xhtml");
		ofs << "<html xmlns=\"http://www.w3.org/1999/xhtml\"><body>";
		g.osHtml = &ofs;

		g.addEdge(v[1], v[1], IMR(1, 1, 1, 1, IV_(3,  1, -1,  0)));
		g.addEdge(v[2], v[2], IMR(1, 1, 1, 1, IV_(3,  1, -1,  0)));

		g.addEdge(v[1], v[2], IMR(1, 1, 1, 1, IV_(3,  1, -1,  0)));
		g.addEdge(v[2], v[1], IMR(1, 1, 1, 1, IV_(3,  1, -1,  0)));

		g.addEdge(v[1], v[0], IMR(1, 3, 1, 2, IV_(4,  1, -1,  0,  0),
				                              IV_(4,  1,  0, -1,  0)));
		g.addEdge(v[0], v[1], IMR(3, 1, 1, 2, IV_(4,  1,  0, -1,  0),
		         							  IV_(4,  0,  1, -1,  0)));
		g.addEdge(v[2], v[0], IMR(1, 3, 1, 2, IV_(4,  1, -1,  0,  0),
				                              IV_(4,  1,  0, -1,  0)));
		g.addEdge(v[0], v[2], IMR(3, 1, 1, 2, IV_(4,  1,  0, -1,  0),
		 		                              IV_(4,  0,  1, -1,  0)));

		g.addEdge(v[0], v[0], IMR(3, 3, 1, 2, IV_(5,  1,  0, -1,  0,  0),
				                              IV_(5,  0,  1,  0, -1,  0)));

		RADGraph* radg = (RADGraph*) g.calcConjunctiveTransitiveClosure();

		vector<RADVertex* > w = vector<RADVertex* >();
		for (unsigned int q = 0; q < v.size(); q++) w.push_back(&radg->getVertex(q));

		vector<vector<vector<RADEdge*> > > resultEdges = extractResultEdges(radg, w);

		ofs << "</body></html>"; ofs.close();
	}*/
}

cute::suite* Test_RADGraph_runSuite(){
	cute::suite &s = *(new cute::suite("RADGraph"));

//	s.push_back(CUTE(Test_RADGraph::CRADGraph_TransitiveClosure));

	return &s;
}
