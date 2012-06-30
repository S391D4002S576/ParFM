#ifndef GRAPH_H_
#define GRAPH_H_

#include <vector>
using std::vector;

#include <string>
using std::string;

#include "../cute/cute.h"

#include <sstream>
#include <stdlib.h>
#include <fstream>
#include <iostream>
using std::ostream;

#include "../basicmath/scalar/Integer.h"
#include "../core/util/DebugStream.h"

#include "../Declarations.h"
#define VXI(VariableName, Graf) ITT(vector<Graf::Vertex* >, VariableName, Graf.getVertices())

namespace AlgoTrans {

	template <class VD, class ED> class CVertex;
	template <class VD, class ED> class CEdge;

	template <class VD, class ED> class CGraph {
		public:
			typedef CVertex<VD, ED> Vertex;
			typedef CEdge<VD, ED> Edge;
			typedef CGraph<VD, ED> Graph;
		protected:
			vector<Vertex*> vertices;
			vector<Edge*> edges;


		public:
			ostream* osLatex;
			DebugOutStream* osHtml;

			vector<Vertex* >& getVertices() { return vertices; }
			const vector<Vertex* >& getVertices() const { return vertices; }

			CGraph() : osLatex(NULL), osHtml(NULL) {};
			CGraph(const CGraph& other);
			virtual ~CGraph();

			virtual Graph* newInstance() const { return new CGraph<VD, ED>(); }

			CGraph& operator = (const CGraph& other);

			const int getVertexCount() const { return vertices.size(); };
			const int getEdgeCount() const { return edges.size(); };

			Vertex& operator () (int index) { return getVertex(index); }
			const Vertex& operator () (int index) const { return getVertex(index); }
			Vertex& getVertex(int index) { return *vertices[index]; };
			const Vertex& getVertex(int index) const { return *vertices[index]; };
			Edge& getEdge(int index) { return *edges[index]; };

			int registerVertex(Vertex& v) {
				vertices.push_back(&v);

				return getVertexCount() - 1;
			};

			typedef typename vector<Vertex*>::iterator iterator;
			iterator begin() { return getVertices().begin(); };
			iterator end() { return getVertices().end(); };

			Vertex& addVertex();
			Vertex& addVertex(const VD& data);
			Vertex& addClonedVertex(const Vertex& originalVertex) { return addVertex(originalVertex.getData()); }

			Edge& addEdge(Vertex& fromVertex, Vertex& toVertex);
			Edge& addEdge(Vertex& fromVertex, Vertex& toVertex, ED* data); // Transfers ownership of data to the edge
			Edge& addEdge(Vertex& fromVertex, Vertex& toVertex, const ED& data); // Takes a copy of the data
			Edge& addEdge(Vertex* fromVertex, Vertex* toVertex, const ED& data); // Takes a copy of the data

			void duplicateVertices(CGraph<VD, ED>& otherG) { transferVertices(otherG); }
			template <class OVD, class OED> void transferVertices(CGraph<OVD, OED>& otherG);

			int registerEdge(Edge* e);

			const vector<Edge*> getEdgesBetween(Vertex& fromVertex, Vertex& toVertex);
			const vector<Edge*> getEdgesBetween(const Vertex& fromVertex, const Vertex& toVertex) const;
			const vector<Edge*> getEdgesBetween(int fromVertexIndex, int toVertexIndex);
			Edge& getUniqueEdgeBetween(Vertex& fromVertex, Vertex& toVertex);
			Edge& getUniqueEdgeBetween(int fromVertexIndex, int toVertexIndex);


			template <class Operations> Graph* floydWarshall(Operations& operations) const;
			template <class Operations> Graph* combineConcurrentEdges(Operations& operations);

			Graph* calcTransitiveClosure();
			Graph* calcUnigraphClosure();

			Graph* calcConjunctiveTransitiveClosure();
			Graph* calcJoinTransitiveClosure();
			Graph* joinConcurrentEdges();

			Graph* calcApproximateDisjunctiveCombinationReduction();

			vector<vector<Vertex*> > getConnectedComponents();
			vector<vector<VD> > getConnectedComponentsData();
			vector<Graph*> splitIntoConnectedComponents();
			vector<Graph> getVertexSetInducedSubGraphs(vector<vector<Vertex*> > comps);
			CGraph<Graph, bool> getVertexSetGraphInducedGraph(CGraph<vector<Vertex*>, bool > comps);

			CGraph<Graph, bool> getStronglyConnectedDecomposition();
			CGraph<vector<CVertex<VD, ED>* >, bool> getStronglyConnectedDecompositionSCCs();
			CGraph<vector<VD>, bool> getStronglyConnectedDecompositionSCCsData();

			Graph subGraph(vector<Vertex*> vertices);

			std::string getDotCode();

			void saveToDotStream(ostream& os);
			void saveToSVGStream(ostream& os);
			std::string saveToSVGString(int ix);
			void saveToDotFile(const std::string& fileName);
			void saveToSVG(const std::string& fileName);
	};

	enum UGraph_Direction {preceding, succeeding, duplex};

	template <class VD, class ED> class CVertex {
		typedef CVertex<VD, ED> Vertex;
		typedef CEdge<VD, ED> Edge;
		typedef CGraph<VD, ED> Graph;
	protected:
		std::string name;

		Graph& graph;
		int index;
		VD* data;

		vector<Edge*> outgoingEdges;
		vector<Edge*> incomingEdges;
	public:
		CVertex(Graph& iGraph);
		CVertex(Graph& iGraph, const VD& data);

		~CVertex() { delete data; }

		const vector<Edge*>& getOutgoingEdges() const { return outgoingEdges; };
		const vector<Edge*>& getIncomingEdges() const { return incomingEdges; };

		vector<Edge*>& getRelatedEdges(UGraph_Direction direction) {
			switch (direction) {
				case preceding: return getIncomingEdges();
				case succeeding: return getOutgoingEdges();
				//default: throw Exception(); XXX
			};
		};

		Edge* getIncomingEdge(const Vertex& predecessorVertex);

		void addRelatedEdge(UGraph_Direction direction, Edge& edge) { getRelatedEdges(direction).push_back(&edge); }
		void addOutgoingEdge(Edge& edge) { addRelatedEdge(succeeding, edge); }
		void addIncomingEdge(Edge& edge) { addRelatedEdge(preceding, edge); }

		int getPredecessorCount() { return getPredecessors().size(); }
		int getSuccessorCount() { return getSuccessors().size(); }

		vector<Vertex*> getPredecessors() { return getRelatedVertices(preceding); };
		vector<Vertex*> getSuccessors() { return getRelatedVertices(succeeding); };

		vector<Vertex*> getRelatedVertices(UGraph_Direction graph_direction);

		int getIndex() const { return index; };

		void registerIncomingEdge(Edge* e);
		void registerOutgoingEdge(Edge* e);

		VD& getData() { return *data; }
		const VD& getData() const { return *data; }

		void setData(const VD newData) {
			if (data == NULL) {
				data = new VD(newData);
			} else {
				*data = newData;
			}
		}

		void setData(VD* newData) {
			if (data != NULL) delete data;

			data = newData;
		}

		const Graph& getGraph() const { return graph; }

		void setName(std::string newName) { name = newName; }
		std::string getName() const { return name; }
		std::string getNonEmptyName() const { return (getName() != "") ? getName() : (CInteger(getIndex()).toString()); }
	};

	template <class VD, class ED> class CEdge {
		template <class VC, class EC> friend class CGraph;

		typedef CVertex<VD, ED> Vertex;
		typedef CEdge<VD, ED> Edge;
		typedef CGraph<VD, ED> Graph;
	protected:
		ED* data;

		Graph& graph;
		Vertex& fromVertex;
		Vertex& toVertex;
	public:
		CEdge(Graph& iGraph, Vertex& iFromVertex, Vertex& iToVertex);
		CEdge(Graph& iGraph, Vertex& iFromVertex, Vertex& iToVertex, const ED& data);
		CEdge(Graph& iGraph, Vertex& iFromVertex, Vertex& iToVertex, ED* data);

		~CEdge() { delete data; }

		/** Data functions **/
		Vertex& getFromVertex() { return fromVertex; };
		Vertex& getToVertex() { return toVertex; };

		ED& getData() { return *data; }
		const ED& getData() const { return *data; }
	};
}

#include "Graph.cpp"
#include "Vertex.cpp"
#include "Edge.cpp"

namespace AlgoTrans {
	template <class VD, class ED>
	vector<CGraph<VD, ED> > CGraph<VD, ED>::getVertexSetInducedSubGraphs(vector<vector<Vertex*> > comps) {
		vector<CGraph<VD, ED> > result = vector<CGraph<VD, ED> >();

		ITT(vector<vector<Vertex*> >, c, comps) result.push_back(subGraph(*c));

		return result;
	}

	template <class VD, class ED>
	CGraph<CGraph<VD, ED>, bool> CGraph<VD, ED>::getVertexSetGraphInducedGraph(CGraph<vector<Vertex*>, bool> comps) {
		typedef vector<CVertex<vector<Vertex*>, bool>* > SCCVertexV;
		typedef CVertex<CGraph<VD, ED>, bool> SCCVertexG;
		typedef vector<CEdge<vector<Vertex* >, bool>* > SCCEdgeV;

		CGraph<CGraph<VD, ED>, bool> result = CGraph<CGraph<VD, ED>, bool>();

		ITT(SCCVertexV, scc, comps.getVertices()) result.addVertex(subGraph((*scc)->getData()));

		//SCCVertexG sccGA = result.getVertices().begin();
		ITT(SCCVertexV, sccA, comps.getVertices()) {
			//SCCVertexG sccGB = result.getVertices().begin();
			ITT(SCCVertexV, sccB, comps.getVertices()) {
				vector<CEdge<vector<Vertex* >, bool>* > edges = comps.getEdgesBetween((*sccA)->getIndex(), (*sccB)->getIndex());

				ITT(SCCEdgeV, e, edges) result.addEdge(result((*sccA)->getIndex()), result((*sccB)->getIndex()), (*e)->getData());

				//++sccB;
			}
			//++sccA;
		}

		return result;
	}

	template <class VD, class ED>
	CGraph<VD, ED> CGraph<VD, ED>::subGraph(vector<Vertex*> vertices) {
		typedef typename vector<Vertex*>::const_iterator VertexIterator;
		typedef typename vector<Edge*>::const_iterator EdgeIterator;

		CGraph<VD, ED> result = CGraph<VD, ED>();
		ITTT(VertexIterator, v, vertices) {
			//(*v)->getData().print();
			//CInteger((*v)->getIndex()).print();
			result.addClonedVertex(**v);
			//result.getVertex(result.getVertexCount() - 1).getData().print();
		}

		VertexIterator vR = result.getVertices().begin();
		ITTT(VertexIterator, v, vertices) {
			VertexIterator wR = result.getVertices().begin();
			ITTT(VertexIterator, w, vertices) {
				vector<Edge*> edges = getEdgesBetween(**v, **w);
				ITTT(EdgeIterator, e, edges) result.addEdge(**vR, **wR, (*e)->getData());
				++wR;
			}
			++vR;
		}

		result.osHtml = osHtml;

		return result;
	}

	template <class VD, class ED> template <class OVD, class OED>
	void CGraph<VD, ED>::transferVertices(CGraph<OVD, OED>& otherG) {
		ITT(vector<Vertex* >, v, otherG.getVertices()) addVertex((*v)->getData());
	}

	template <class VD, class ED>
	std::string CGraph<VD, ED>::getDotCode() {
		std::string result = "digraph GraphName {\n";

		// Nodes
		result += "node [shape=circle];";
		ITT(vector<Vertex*>, v, getVertices()) result += " " + (*v)->getNonEmptyName() + ";";
		result += "\n";

		// Edges
		ITT(vector<Vertex*>, v, getVertices()) {
			const vector<Edge*>& outgEdges = (*v)->getOutgoingEdges();
			for (unsigned int q = 0; q < outgEdges.size(); q++) { Edge* e = outgEdges[q];
				result += (e->getFromVertex().getNonEmptyName()) + "->" + (e->getToVertex().getNonEmptyName()) + ";\n";
			}
		}

		result += "overlap = false\n";

		return result + "}";
	}

	template <class VD, class ED>
	void CGraph<VD, ED>::saveToDotStream(ostream& os) {
		os << getDotCode();
	}

	template <class VD, class ED>
	void CGraph<VD, ED>::saveToDotFile(const std::string& fileName) {
		std::ofstream ofs(fileName.c_str());

		saveToDotStream(ofs);

		ofs.close();
	}

	template <class VD, class ED>
	void CGraph<VD, ED>::saveToSVG(const std::string& fileName) {
		string tmpDotCode = fileName + ".dot";
		saveToDotFile(tmpDotCode);

		string command = "neato -Tsvg " + tmpDotCode + " > " + fileName;
		system(command.c_str());
		remove(tmpDotCode.c_str());
	}

	template <class VD, class ED>
	void CGraph<VD, ED>::saveToSVGStream(ostream& os) { os << saveToSVGString(); }

	template <class VD, class ED>
	std::string CGraph<VD, ED>::saveToSVGString(int ix) {
		std::string tmpFile = "tempSVGFile" + CInteger(ix).toString() + ".svg"; // not thread safe -- not even remotely decent --> find tempFileName-getting function
		saveToSVG("examples/" + tmpFile);

		std::ifstream ifs(tmpFile.c_str(), std::ios::binary);
		ifs.seekg(0, std::ios::end);
		int length = ifs.tellg();
		ifs.seekg(0, std::ios::beg);
		char* buffer = new char[length + 2];
		buffer[length] = 0; buffer[length - 1] = 0;
		ifs.read(buffer, length);
		std::string svgString = buffer;
		delete[] buffer;

		unsigned int svgStartPos = svgString.find("<svg");
		//printf("%i -- %i\n", svgStartPos, svgString.npos);
		//if (svgStartPos != svgString.npos) svgString = svgString.substr(svgStartPos);
		ifs.close();

		//remove(tmpFile.c_str());
		return "<img src=\"" + tmpFile + "\">";
	}
}

#endif /*GRAPH_H_*/
