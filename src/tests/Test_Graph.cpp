#include "../cute/cute.h"
#include "../cute/ide_listener.h"
#include "../cute/cute_suite.h"

#include "../partyql/graph/Graph.h"
#include "../partyql/StandardOperations.h"
#include "../partyql/basicmath/scalar/Integer.h"

#include <fstream>
using std::ofstream;

using namespace AlgoTrans;

namespace Test_Graph {

typedef CVertex<void*, CInteger> Vertex;
typedef CEdge<void*, CInteger> Edge;
typedef CGraph<void*, CInteger> Graph;

vector<vector<vector<Edge*> > > extractResultEdges(CGraph<void*, CInteger>* gFW, vector<Vertex* >& w) {
	vector<vector<vector<Edge*> > > resultEdges
		= vector<vector<vector<Edge*> > >(w.size(), vector<vector<Edge*> >(w.size()));

	for (unsigned int q = 0; q < w.size(); q++) {
		for (unsigned int r = 0; r < w.size(); r++) {
			resultEdges[q][r] = gFW->getEdgesBetween(*w[q], *w[r]);
		}
	}

	return resultEdges;
}

CShortestPathOperations<CInteger>* shortestPathOps;

void CGraph_InitOperations() {
	shortestPathOps = new CShortestPathOperations<CInteger>();
}

void CGraph_FloydWarshall_GetEdgesBetween() {
	CGraph_InitOperations();

	Graph g = Graph();
	vector<Vertex* > v = vector<Vertex* >(4);
	for (unsigned int q = 0; q < v.size(); q++) v[q] = &g.addVertex();

	Edge& e01 = g.addEdge(*v[0], *v[1], CInteger(1));
	/* Edge& e12a = */ g.addEdge(*v[1], *v[2], CInteger(1));
	/* Edge& e12b = */ g.addEdge(*v[1], *v[2], CInteger(2));
	Edge& e23 = g.addEdge(*v[2], *v[3], CInteger(1));

	vector<vector<vector<Edge*> > > resultEdges = extractResultEdges(&g, v);

	ASSERT((resultEdges[0][1].size() == 1) && (resultEdges[0][1][0] == &e01));
	ASSERT((resultEdges[1][0].size() == 0));
	ASSERT((resultEdges[1][2].size() == 2));
	ASSERT((resultEdges[2][1].size() == 0));
	ASSERT((resultEdges[0][2].size() == 0));
	ASSERT((resultEdges[0][3].size() == 0));
	ASSERT((resultEdges[2][3].size() == 1) && (resultEdges[2][3][0] == &e23));
}

/*void CGraph_FloydWarshall_ShortestPath_1() {
	CGraph_InitOperations();

	Graph g = Graph();
	vector<Vertex* > v = vector<Vertex* >(4);
	for (unsigned int q = 0; q < v.size(); q++) v[q] = &g.addVertex();

	g.addEdge(*v[0], *v[1], CInteger(1));
	g.addEdge(*v[1], *v[2], CInteger(1));
	g.addEdge(*v[2], *v[3], CInteger(1));

	ofstream ofs("examples/CGraph_FloydWarshall_ShortestPath_1.xhtml");
	ofs << "<html xmlns=\"http://www.w3.org/1999/xhtml\"><body>";
	g.osHtml = &ofs;

	Graph* gFW = g.floydWarshall<CShortestPathOperations<CInteger> >(*shortestPathOps);

	ofs << "</body></html>"; ofs.close();

	vector<Vertex* > w = vector<Vertex* >();
	for (unsigned int q = 0; q < v.size(); q++) w.push_back(&gFW->getVertex(q));

	vector<vector<vector<Edge*> > > resultEdges = extractResultEdges(gFW, w);

	vector<Edge*> ve = resultEdges[0][1];
	ASSERT((resultEdges[0][1].size() == 1) && (resultEdges[0][1][0]->getData() == CInteger(1)));
	ASSERT((resultEdges[1][0].size() == 0));
	ASSERT((resultEdges[1][2].size() == 1) && (resultEdges[1][2][0]->getData() == CInteger(1)));
	ASSERT((resultEdges[2][1].size() == 0));
	ASSERT((resultEdges[0][2].size() == 1) && (resultEdges[0][2][0]->getData() == CInteger(2)));
	ASSERT((resultEdges[0][3].size() == 1) && (resultEdges[0][3][0]->getData() == CInteger(3)));
}

void CGraph_FloydWarshall_ShortestPath_2() {
	CGraph_InitOperations();

	CGraph<void*, CInteger> g = CGraph<void*, CInteger>();
	vector<Vertex* > v = vector<Vertex* >(4);
	for (unsigned int q = 0; q < v.size(); q++) v[q] = &g.addVertex();

	g.addEdge(*v[0], *v[1], CInteger(1));
	g.addEdge(*v[1], *v[2], CInteger(4));
	g.addEdge(*v[2], *v[3], CInteger(1));
	g.addEdge(*v[0], *v[2], CInteger(3));

	Graph* gFW = g.floydWarshall<CShortestPathOperations<CInteger> >(*shortestPathOps);

	vector<Vertex* > w = vector<Vertex* >();
	for (unsigned int q = 0; q < v.size(); q++) w.push_back(&gFW->getVertex(q));

	vector<vector<vector<Edge*> > > resultEdges = extractResultEdges(gFW, w);

	ASSERT((resultEdges[0][1].size() == 1) && (resultEdges[0][1][0]->getData() == CInteger(1)));
	ASSERT((resultEdges[1][0].size() == 0));
	ASSERT((resultEdges[1][2].size() == 1) && (resultEdges[1][2][0]->getData() == CInteger(4)));
	ASSERT((resultEdges[2][1].size() == 0));
	ASSERT((resultEdges[0][2].size() == 1) && (resultEdges[0][2][0]->getData() == CInteger(3)));
	ASSERT((resultEdges[0][3].size() == 1) && (resultEdges[0][3][0]->getData() == CInteger(4)));
}

void CGraph_FloydWarshall_ShortestPath_3() {
	CGraph_InitOperations();

	CGraph<void*, CInteger> g = CGraph<void*, CInteger>();
	vector<Vertex* > v = vector<Vertex* >(4);
	for (unsigned int q = 0; q < v.size(); q++) v[q] = &g.addVertex();

	g.addEdge(*v[0], *v[1], CInteger(1));
	g.addEdge(*v[1], *v[2], CInteger(2));
	g.addEdge(*v[2], *v[3], CInteger(1));
	g.addEdge(*v[0], *v[2], CInteger(5));

	Graph* gFW = g.floydWarshall<CShortestPathOperations<CInteger> >(*shortestPathOps);

	vector<Vertex* > w = vector<Vertex* >();
	for (unsigned int q = 0; q < v.size(); q++) w.push_back(&gFW->getVertex(q));

	vector<vector<vector<Edge*> > > resultEdges = extractResultEdges(gFW, w);

	ASSERT((resultEdges[0][1].size() == 1) && (resultEdges[0][1][0]->getData() == CInteger(1)));
	ASSERT((resultEdges[1][0].size() == 0));
	ASSERT((resultEdges[1][2].size() == 1) && (resultEdges[1][2][0]->getData() == CInteger(2)));
	ASSERT((resultEdges[2][1].size() == 0));
	ASSERT((resultEdges[0][2].size() == 1) && (resultEdges[0][2][0]->getData() == CInteger(3)));
	ASSERT((resultEdges[0][3].size() == 1) && (resultEdges[0][3][0]->getData() == CInteger(4)));
}*/

}

cute::suite* Test_Graph_runSuite(){
	cute::suite& s = *(new cute::suite("Graph"));

	s.push_back(CUTE(Test_Graph::CGraph_FloydWarshall_GetEdgesBetween));
	/*s.push_back(CUTE(Test_Graph::CGraph_FloydWarshall_ShortestPath_1));
	s.push_back(CUTE(Test_Graph::CGraph_FloydWarshall_ShortestPath_2));
	s.push_back(CUTE(Test_Graph::CGraph_FloydWarshall_ShortestPath_3));*/

	return &s;
}

