#ifndef GRAPH_CPP_
#define GRAPH_CPP_

#include <string>

#include "Graph.h"
#include "../Structures.h"

namespace AlgoTrans {
	template <class VD, class ED>
	int CGraph<VD, ED>::registerEdge(CEdge<VD, ED>* e) {
		edges.push_back(e);

		e->getToVertex().registerIncomingEdge(e);
		e->getFromVertex().registerOutgoingEdge(e);

		return getEdgeCount() - 1;
	}

	template <class VD, class ED> CVertex<VD, ED>& CGraph<VD, ED>::addVertex() {
		return *(new CVertex<VD, ED>(*this));
	}

	template <class VD, class ED> CVertex<VD, ED>& CGraph<VD, ED>::addVertex(const VD& data) {
		return *(new CVertex<VD, ED>(*this, data));
	}

	template <class VD, class ED>
	CEdge<VD, ED>& CGraph<VD, ED>::addEdge(CVertex<VD, ED>& fromVertex, CVertex<VD, ED>& toVertex) {
		return *(new CEdge<VD, ED>(*this, fromVertex, toVertex));
	}

	template <class VD, class ED> // Transfers ownership of the data
	CEdge<VD, ED>& CGraph<VD, ED>::addEdge(CVertex<VD, ED>& fromVertex, CVertex<VD, ED>& toVertex, ED* data) {
		CEdge<VD, ED>&  edge = addEdge(fromVertex, toVertex);

		edge.data = data;

		return edge;
	}

	template <class VD, class ED> // Takes a copy of the data
	CEdge<VD, ED>& CGraph<VD, ED>::addEdge(CVertex<VD, ED>& fromVertex, CVertex<VD, ED>& toVertex, const ED& data) {
		return *(new CEdge<VD, ED>(*this, fromVertex, toVertex, data));
	}

	template <class VD, class ED> // Takes a copy of the data
	CEdge<VD, ED>& CGraph<VD, ED>::addEdge(CVertex<VD, ED>* fromVertex, CVertex<VD, ED>* toVertex, const ED& data) {
		return *(new CEdge<VD, ED>(*this, *fromVertex, *toVertex, data));
	}

	template <class VD, class ED> const vector<CEdge<VD, ED>*> CGraph<VD, ED>::getEdgesBetween(CVertex<VD, ED>& fromVertex, CVertex<VD, ED>& toVertex) {
		vector<CEdge<VD, ED>*> result = vector<CEdge<VD, ED>*>();

		const vector<CEdge<VD, ED>*>& outEdges = fromVertex.getOutgoingEdges();
		for (unsigned int e = 0; e < outEdges.size(); e++) {
			if (&outEdges[e]->getToVertex() == &toVertex) result.push_back(outEdges[e]);
		}

		return result;
	}

	template <class VD, class ED> const vector<CEdge<VD, ED>*> CGraph<VD, ED>::getEdgesBetween(const CVertex<VD, ED>& fromVertex, const CVertex<VD, ED>& toVertex) const {
		vector<CEdge<VD, ED>*> result = vector<CEdge<VD, ED>*>();

		const vector<CEdge<VD, ED>*>& outEdges = fromVertex.getOutgoingEdges();
		for (unsigned int e = 0; e < outEdges.size(); e++) {
			if (&outEdges[e]->getToVertex() == &toVertex) result.push_back(outEdges[e]);
		}

		return result;
	}

	template <class VD, class ED> const vector<CEdge<VD, ED>*> CGraph<VD, ED>::getEdgesBetween(int fromVertexIndex, int toVertexIndex) {
		return getEdgesBetween(getVertex(fromVertexIndex), getVertex(toVertexIndex));
	}

	template <class VD, class ED>
	CGraph<VD, ED>* CGraph<VD, ED>::calcConjunctiveTransitiveClosure() {
		typename ED::ConjunctiveOperations ops = typename ED::ConjunctiveOperations();

		return floydWarshall<typename ED::ConjunctiveOperations>(ops);
	}

	template <class VD, class ED>
	CGraph<VD, ED>* CGraph<VD, ED>::calcUnigraphClosure() {
		typedef typename ED::TransitiveClosureOperations Ops;

		Ops ops = Ops();

		return combineConcurrentEdges<Ops>(ops);
	}

	template <class VD, class ED>
	CGraph<VD, ED>* CGraph<VD, ED>::joinConcurrentEdges() {
		typedef typename ED::DisjunctiveHullOperations Operations;

		Operations ops = Operations();

		return combineConcurrentEdges<Operations>(ops);
	}

	template <class VD, class ED>
	CGraph<VD, ED>* CGraph<VD, ED>::calcJoinTransitiveClosure() {
		typename ED::DisjunctiveHullOperations ops = typename ED::DisjunctiveHullOperations();

		return floydWarshall<typename ED::DisjunctiveHullOperations>(ops);
	}

	template <class VD, class ED>
	CGraph<VD, ED>* CGraph<VD, ED>::calcApproximateDisjunctiveCombinationReduction() {
		typename ED::DisjunctiveHullOperations ops = typename ED::DisjunctiveHullOperations();

		return combineConcurrentEdges<typename ED::DisjunctiveHullOperations>(ops);
	}

	template <class VD, class ED>
	template <class Operations>
	CGraph<VD, ED>* CGraph<VD, ED>::combineConcurrentEdges(Operations& operations) {
		CGraph<VD, ED>& resultG = *(newInstance());
		resultG.osHtml = osHtml;
		if (osHtml) *osHtml << "<b>Combining concurrent edges...</b></br>";

		for (int p = 0; p < getVertexCount(); p++) resultG.addVertex(getVertex(p).getData());

		// Generate a matrix for easy access to the edge data
		vector<vector<ED*> > edgeData = vector<vector<ED*> >(getVertexCount(), vector<ED*>(getVertexCount(), NULL));
		for (int p = 0; p < getVertexCount(); p++) { const CVertex<VD, ED>& vxP = getVertex(p);
			const vector<CEdge<VD, ED>*>& outgEdges = vxP.getOutgoingEdges();
			for (int q = 0; q < getVertexCount(); q++) { const CVertex<VD, ED>& vxQ = getVertex(q);
				for (unsigned int e = 0; e < outgEdges.size(); e++) { CEdge<VD, ED>& egE = *outgEdges[e];
					if (&egE.getToVertex() == &vxQ) {
						edgeData[p][q] = new ED((edgeData[p][q] != NULL)
								       ? Operations::calcCombination(*edgeData[p][q], egE.getData())
									   : egE.getData());
					}
				}
			}
		}

		// Dump resulting edge data in graph
		for (int p = 0; p < getVertexCount(); p++) { CVertex<VD, ED>& vxP = resultG.getVertex(p);
			for (int q = 0; q < getVertexCount(); q++) { CVertex<VD, ED>& vxQ = resultG.getVertex(q);
				if (edgeData[p][q] != NULL) resultG.addEdge(vxP, vxQ, edgeData[p][q]);
			}
		}

		return &resultG;
	}

	template <class VD, class ED>
	CGraph<VD, ED>::CGraph(const CGraph<VD, ED>& other) {
		if (this != &other) {
			osHtml = other.osHtml;
			for (int q = 0; q < other.getVertexCount(); q++) {
				addVertex(other.getVertex(q).getData());
			}
			for (int q = 0; q < other.getVertexCount(); q++) {
				for (int r = 0; r < other.getVertexCount(); r++) {
					const vector<Edge*> edges = other.getEdgesBetween(other.getVertex(q), other.getVertex(r));
					for (unsigned int z = 0; z < edges.size(); z++) { Edge* e = edges[z];
						addEdge(getVertex(q), getVertex(r), e->getData());
					}
				}
			}
		}
	}

	template <class VD, class ED>
	CGraph<VD, ED>::~CGraph() {
		ITT(vector<Vertex*>, v, vertices) delete *v;
		ITT(vector<Edge*>, e, edges) delete *e;
	};

	template <class VD, class ED>
	template <class Operations>
	CGraph<VD, ED>* CGraph<VD, ED>::floydWarshall(Operations& operations) const {
		CGraph<VD, ED>& resultG = *(newInstance());

		if (osHtml != NULL) *osHtml << "<b>Starting Floyd-Warshall:</b><br/>";

		for (int p = 0; p < getVertexCount(); p++) { const CVertex<VD, ED>& vxP = getVertex(p);
			resultG.addVertex(vxP.getData());
		}

		//H("<b>Converting multi-graph to uni-graph:</b>");
		//H("<table><tr><td>LS</td><td>RS</td><td>E</td><td>edgeData[p][q]</td><td>egE.getData()</td><td>Result</td></tr>");
		// Generate a matrix for easy access to the edge data
		// We already combine parallel edges here for multigraphs
		vector<vector<ED*> > edgeData = vector<vector<ED*> >(getVertexCount(), vector<ED*>(getVertexCount(), NULL));
		vector<vector<ED*> > edgeStartData = vector<vector<ED*> >(getVertexCount(), vector<ED*>(getVertexCount(), NULL));
		vector<vector<ED*> > edgeData2 = vector<vector<ED*> >(getVertexCount(), vector<ED*>(getVertexCount(), NULL));
		for (int p = 0; p < getVertexCount(); p++) { const CVertex<VD, ED>& vxP = getVertex(p);
			for (int q = 0; q < getVertexCount(); q++) { const CVertex<VD, ED>& vxQ = getVertex(q);
				const vector<CEdge<VD, ED>*>& edges = getEdgesBetween(vxP, vxQ);
			    if (edges.size() > 0) edgeData[p][q] = new ED(edges[0]->getData());
			    /*
				for (unsigned int e = 0; e < outgEdges.size(); e++) { CEdge<VD, ED>& egE = *outgEdges[e];
					if (&egE.getToVertex() == &vxQ) {
						H("<tr><td>"); H(CInteger(p).toString()); H("</td><td>"); H(CInteger(q).toString()); H("</td><td>"); H(CInteger(e).toString());
						H("</td><td>"); H(edgeData[p][q]->toStringHtml()); H("</td><td>"); H(egE.getData().toStringHtml()); H("</td>");

						edgeData[p][q] = new ED((edgeData[p][q] != NULL) ? Operations::calcCombination(*edgeData[p][q], egE.getData())
														 : egE.getData());
						edgeStartData[p][q] = edgeData[p][q];

						H("<td>"); H(edgeStartData[p][q]->toStringHtml()); H("</td></tr>");
					}
				}*/
			}
		}
		//if (osHtml != NULL) *osHtml << "</table>";

		int dims = 1;
		vector<string>* dimNamesI = NULL;
		vector<string>* dimNamesIp = NULL;
		vector<string>* dimNamesIIp = NULL;
		dimNamesI = new vector<string>();
		dimNamesI->push_back("p");
		if (dims == 2) dimNamesI->push_back("s");

		dimNamesIp = new vector<string>();
		dimNamesIp->push_back("q");
		if (dims == 2) dimNamesIp->push_back("t");

		dimNamesIIp = new vector<string>();
		dimNamesIIp->push_back("p");
		if (dims == 2) dimNamesIIp->push_back("s");
		dimNamesIIp->push_back("q");
		if (dims == 2) dimNamesIIp->push_back("t");

		if (osHtml != NULL) *osHtml << "<b>After conversion from multi-graph to unigraph:</b><table>";
		if (osHtml != NULL) *osHtml << "<tr><td></td>";
		for (int q = 0; q < getVertexCount(); q++) {
			H("<td>Statement "); H(CInteger(q).toString()); H("</td>");
		}
		if (osHtml != NULL) *osHtml << "</tr>";
		for (int p = 0; p < getVertexCount(); p++) {
			H("<tr><td>Statement "); H(CInteger(p).toString()); H("</td>");
			for (int q = 0; q < getVertexCount(); q++) {
				H("<td><table><tr><td><center>");
				H(edgeData[p][q]->toStringHtml());
				H("</center></td><td></td></tr></table></td>");
			}
			if (osHtml != NULL) *osHtml << "</tr>";
		}
		if (osHtml != NULL) *osHtml << "</table>";

		// Actual FW
		for (int k = 0; k < getVertexCount(); k++) {
			if (edgeData[k][k] != NULL) { *edgeData[k][k] = edgeData[k][k]->getTransitiveClosure(); }
			for (int l = 0; l < getVertexCount(); l++) {
				for (int m = 0; m < getVertexCount(); m++) {
					if (edgeData2[l][m] != NULL) delete edgeData2[l][m];
					edgeData2[l][m] = (edgeData[l][m] == NULL) ? NULL : new ED(*edgeData[l][m]);
				}
			}

			H("<b>FW-Core processing vertex "); H(CInteger(k).toString()); H(":</b><table>");
			H("<tr><td>Li</td><td>Ri</td><td>L><td>M</td><td>M*</td><td>R</td><td>L o M</td><td>L o M o R</td><td>Original edge</td><td>Combination</td></tr>");
			for (int l = 0; l < getVertexCount(); l++) { if (edgeData[l][k] != NULL) {
				for (int m = 0; m < getVertexCount(); m++) { if (edgeData[k][m] != NULL) {
					ED TC = edgeData[k][k]->getTransitiveClosure();
					ED* concurrD = new ED(Operations::calcTransition(*edgeData[l][k], TC));
					ED* concurrE = new ED(Operations::calcTransition(*concurrD, *edgeData[k][m]));

					H("<tr><td><center>"); H(CInteger(l).toString()); H("</center></td><td><center>");
					H(CInteger(m).toString()); H("</center></td><td><center>");
					H(edgeData[l][k]->toStringHtml()); H("</center></td><td><center>");
					H(edgeData[k][k]->toStringHtml()); H("</center></td><td><center>");
					H(TC.toStringHtml()); H("</center></td><td><center>");
					H(edgeData[k][m]->toStringHtml()); H("</center></td><td><center>");
					H(concurrD->toStringHtml()); H("</center></td><td><center>");
					H(concurrE->toStringHtml()); H("</center></td><td><center>");
					H((edgeData[l][m] != NULL) ? edgeData[l][m]->toStringHtml() : std::string("Empty edge")); H("</center></td><td><center>");
					H(((edgeData[l][m] != NULL) ? Operations::calcCombination(*concurrE, *edgeData[l][m]).toStringHtml() : concurrD->toStringHtml()));
					H("</center></td></tr>");

					if (edgeData2[l][m] != NULL) delete edgeData2[l][m];
					edgeData2[l][m] = new ED((edgeData[l][m] != NULL) ? Operations::calcCombination(*concurrE, *edgeData[l][m])
													 : *concurrE);
					delete concurrD;
				} }
			} }
			if (osHtml != NULL) *osHtml << "</table>";

			//ASSERTHYPO(hypo1(k, edgeStartData, edgeData2, edgeValid2));

			for (int l = 0; l < getVertexCount(); l++) {
				for (int m = 0; m < getVertexCount(); m++) {
					if (edgeData[l][m] != NULL) delete edgeData[l][m];
					edgeData[l][m] = (edgeData2[l][m] != NULL) ? new ED(*edgeData2[l][m]) : NULL;
				}
			}

			H("<b>After FW-core step "); H(CInteger(k).toString()); H(":</b>");
			H("<table>");
			H("<tr><td></td>");
			for (int q = 0; q < getVertexCount(); q++) {
				H("<td>Statement "); H(CInteger(q).toString()); H("</td>");
			}
			H("</tr>");
			for (int l = 0; l < getVertexCount(); l++) {
				H("<tr><td>Statement "); H(CInteger(l).toString()); H("</td>");
				for (int m = 0; m < getVertexCount(); m++) { if (edgeData[l][m] != NULL) {
					H("<td><table><tr><td><center>");
					H(edgeData[l][m]->toStringHtml());
					H("</center></td><td></td></tr></table></td>");
				} else {
					H("<td></td>");
				} }
				if (osHtml != NULL) *osHtml << "</tr>";
			}
			if (osHtml != NULL) *osHtml << "</table>";
		}

		// Dump resulting edge data in graph
		for (int p = 0; p < getVertexCount(); p++) { CVertex<VD, ED>& vxP = resultG.getVertex(p);
			for (int q = 0; q < getVertexCount(); q++) { CVertex<VD, ED>& vxQ = resultG.getVertex(q);
				if (edgeData[p][q] != NULL) {
					resultG.addEdge(vxP, vxQ, *edgeData[p][q]);
					delete edgeData[p][q];
				}
			}
		}

		return &resultG;
	}

	template <class VD, class ED>
	CGraph<VD, ED>& CGraph<VD, ED>::operator = (const CGraph<VD, ED>& other) {
		if (this != &other) {
			vertices = other.vertices;
			edges = other.edges;
			ASSERT(false); // Not yet implemented
		}

		return *this;
	}


	template <class VD, class ED>
	CEdge<VD, ED>& CGraph<VD, ED>::getUniqueEdgeBetween(CVertex<VD, ED>& fromVertex, CVertex<VD, ED>& toVertex) {
		vector<CEdge<VD, ED>*> edges = getEdgesBetween(fromVertex, toVertex);

		ASSERT((edges.size() == 1));

		return *edges[0];
	}

	template <class VD, class ED>
	CEdge<VD, ED>& CGraph<VD, ED>::getUniqueEdgeBetween(int fromVertexIndex, int toVertexIndex) {
		return getUniqueEdgeBetween(getVertex(fromVertexIndex), getVertex(toVertexIndex));
	}

	/* Connected Component Finder */
	template <class VD, class ED> class CConnectedComponentFinder {
		typedef CVertex<VD, ED> Vertex;
	private:
		CGraph<VD, ED>& graph;

		vector<Vertex* > vertexList;
		vector<bool> vertexOnStack, vertexVisited;

		vector<Vertex* >* vertexStack;

		void connectedComponentRecursion(int vix) {
			Vertex& vertex = graph(vix);
			vertexVisited[vix] = true;

			vertexStack->push_back(&vertex);

			vector<Vertex* > successors = vertex.getRelatedVertices(duplex);
			ITT(vector<Vertex* >, succV, successors) { int succVix = (*succV)->getIndex();
				if (!vertexVisited[succVix]) connectedComponentRecursion(succVix);
			}
		}
	public:
		vector<vector<Vertex* > > findConnectedComponents() {
			vector<vector<Vertex* > > result;

			vertexVisited = vector<bool>(graph.getVertexCount(), false);

			int nextVertex = graph.getVertexCount();
			while (--nextVertex >= 0) if (!vertexVisited[nextVertex]) {
				result.push_back(vector<Vertex* >());
				vertexStack = &result[result.size() - 1];

				connectedComponentRecursion(nextVertex);
			}

			return result;
		}

		CConnectedComponentFinder(CGraph<VD, ED>& iGraph) : graph(iGraph) { };
	};

	template <class VD, class ED>
	vector<vector<CVertex<VD, ED>* > > CGraph<VD, ED>::getConnectedComponents() {
		return CConnectedComponentFinder<VD, ED>(*this).findConnectedComponents();
	}

	/*template <class VD, class ED>
	vector<vector<int> > CGraph<VD, ED>::getConnectedComponentsByIndex() {
		vector<vector<CVertex<VD, ED>* > > v = getConnectedComponents();

		vector<vector<int> > result = vector<vector<CInteger> >();
		for (unsigned int q = 0; q < v.size(); q++) {
			vector<int> v = vector<int>();
			for (unsigned int f = 0; f < v[q].size(); f++) {
				v.push_back(v[q][f]->getIndex());
			}
			result.push_back(v);
		}

		return v;
	}*/

	template <class VD, class ED>
	vector<vector<VD> > CGraph<VD, ED>::getConnectedComponentsData() {
		vector<vector<CVertex<VD, ED>* > > comps = getConnectedComponents();
		vector<vector<VD> > result = vector<vector<VD> >();
		for (unsigned int q = 0; q < comps.size(); q++) {
			result.push_back(vector<VD>());
			for (unsigned int r = 0; r < comps[q].size(); r++) {
				result[q].push_back(comps[q][r]->getData());
			}
		}

		return result;
	}

	template <class VD, class ED>
	vector<CGraph<VD, ED>* > CGraph<VD, ED>::splitIntoConnectedComponents() {
		vector<vector<Vertex* > > comps = getConnectedComponents();

		vector<Graph*> result = vector<Graph*>();

		// Subcomponent graphs first
		for (typename vector<vector<Vertex* > >::iterator comp = comps.begin(); comp != comps.end(); ++comp) {
			Graph* compGraph = new Graph();
			for (unsigned int q = 0; q < (*comp).size(); q++) compGraph.addVertex((*comp)[q]->getData());
			for (unsigned int q = 0; q < (*comp).size(); q++) {
				for (unsigned int r = 0; r < (*comp).size(); r++) {
					vector<Edge* > edges = getEdgesBetween((*comp)[q], (*comp)[r]);
					for (typename vector<Edge*>::iterator e = edges.begin(); e != edges.end(); ++e) {
						compGraph.addEdge(q, r, (*e)->getData());
					}
				}
			}
			result.push_back(compGraph);
		}

		return result;
	}

	/* Strongly Connected Component Finder: Based on (originally inspired by) Tarjan's algorithm */
	template <class VD, class ED> class CStronglyConnectedComponentFinder {
		typedef CVertex<VD, ED> Vertex;
		typedef typename std::vector<Vertex* >::iterator VertexIt;
	private:
		CGraph<VD, ED>& graph;

		vector<int> lowvisIndex;
		vector<bool> vertexVisited;

		int verticesVisited;

		CGraph<vector<Vertex* >, bool> result;
		CVertex<vector<Vertex* >, bool>* currentComp;
		vector<CVertex<vector<Vertex* >, bool>* > activeComps;
		vector<int> activeCompsLowvis;

		void visit(int vix) {
			Vertex& vertex = graph(vix);
			vertexVisited[vix] = true;

			lowvisIndex[vix] = verticesVisited++;
			int visitIndex = lowvisIndex[vix];

			vector<Vertex*> successors = vertex.getSuccessors();
			ITT(vector<Vertex*>, succV, successors) { int succVix = (*succV)->getIndex();
				if (!vertexVisited[succVix]) visit(succVix);
				lowvisIndex[vix] = std::min(lowvisIndex[vix], lowvisIndex[succVix]);
			}

			if (currentComp == NULL) { currentComp = &result.addVertex(); }
			currentComp->getData().push_back(&vertex);

			if (visitIndex == lowvisIndex[vix]) { // First vertex of a strongly connected component -- last vertex of the component in the (returning of) the recursion
				if (activeComps.size() > 0) { // Connect with succeeding subgraphs in reduced graph
					int rem = 0;
					for (unsigned int q = 0; q < activeComps.size(); q++) {
						activeComps[rem] = activeComps[q];
						activeCompsLowvis[rem] = activeCompsLowvis[q];
						if (activeCompsLowvis[q] > visitIndex) {
							result.addEdge(*currentComp, *(activeComps[q]), true);
						} else rem++;
					}
					activeComps.resize(rem); activeCompsLowvis.resize(rem);
				}
				activeComps.push_back(currentComp);
				activeCompsLowvis.push_back(visitIndex);

				currentComp = NULL;
			}
		}
	public:
		CGraph<vector<Vertex* >, bool> findStronglyConnectedComponents() {
			lowvisIndex = vector<int>(graph.getVertexCount(), -1);
			vertexVisited = vector<bool>(graph.getVertexCount(), false);

			int nextVertex = graph.getVertexCount();

			verticesVisited = 0;
			currentComp = NULL;
			while (--nextVertex >= 0) if (!vertexVisited[nextVertex]) visit(nextVertex);

			return result;
		}

		CGraph<CGraph<VD, ED>, bool> getDecomposedGraph() {
			typedef CVertex<vector<Vertex* >, bool> SCCVertexV;
			typedef CVertex<CGraph<VD, ED>, bool> SCCVertexG;
			typedef vector<CEdge<vector<Vertex* >, bool>* > SCCEdgeV;

			CGraph<vector<Vertex* >, bool> sccGraph = findStronglyConnectedComponents();
			CGraph<CGraph<VD, ED>, bool> result = CGraph<CGraph<VD, ED>, bool>();

			ITT(SCCVertexV, scc, sccGraph.getVertices()) sccGraph.addVertex(graph.subGraph(*scc));

			SCCVertexG sccGA = result.getVertices().begin();
			ITT(SCCVertexV, sccA, sccGraph.getVertices()) {
				SCCVertexG sccGB = result.getVertices().begin();
				ITT(SCCVertexV, sccB, sccGraph.getVertices()) {
					SCCEdgeV edges = sccGraph.getEdgesBetween(*sccA, *sccB);

					ITT(SCCEdgeV, e, edges) result.addEdge(*sccGA, *sccGB, (*e)->getData());
					++sccGB;
				}
				++sccGA;
			}

			return result;
		}

		CStronglyConnectedComponentFinder(CGraph<VD, ED>& iGraph) : graph(iGraph) { };
	};

	template <class VD, class ED>
	CGraph<vector<CVertex<VD, ED>* >, bool> CGraph<VD, ED>::getStronglyConnectedDecompositionSCCs() {
		return CStronglyConnectedComponentFinder<VD, ED>(*this).findStronglyConnectedComponents();
	}

	template <class VD, class ED>
	CGraph<vector<VD>, bool> CGraph<VD, ED>::getStronglyConnectedDecompositionSCCsData() {
		typedef CGraph<vector<CVertex<VD, ED>* >, bool> G;
		typedef CEdge<vector<CVertex<VD, ED>* >, bool> Edge;
		typedef vector<CVertex<VD, ED>* > VectorIt;
		CGraph<vector<CVertex<VD, ED>* >, bool> sccDecomp = getStronglyConnectedDecompositionSCCs();

		CGraph<vector<VD>, bool> result = CGraph<vector<VD>, bool>();

		ITT(G, v, sccDecomp) {
			vector<VD> vect = vector<VD>();
			ITT(VectorIt, m, (*v)->getData()) {
				vect.push_back((*m)->getData());
			}
			result.addVertex(vect);
		}

		ITT(G, v, sccDecomp) {
			ITT(G, w, sccDecomp) {
				vector<Edge*> edges = sccDecomp.getEdgesBetween((*v)->getIndex(), (*w)->getIndex());
				ITT(vector<Edge*>, e, edges) result.addEdge(result((*v)->getIndex()), result((*w)->getIndex()), (*e)->getData());
			}
		}

		return result;
	}

	template <class VD, class ED>
	CGraph<CGraph<VD, ED>, bool> CGraph<VD, ED>::getStronglyConnectedDecomposition() {
		return CStronglyConnectedComponentFinder<VD, ED>(*this).getDecomposedGraph();
	}
}

#endif /* GRAPH_CPP_ */
