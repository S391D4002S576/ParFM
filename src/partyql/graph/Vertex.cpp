#ifndef VERTEX_CPP_
#define VERTEX_CPP_

#include "Graph.h"
#include "../Structures.h"

namespace AlgoTrans{
	template <class VD, class ED> CVertex<VD, ED>::CVertex(CGraph<VD, ED>& iGraph): graph(iGraph) {
		index = graph.registerVertex(*this);

		data = new VD();
	}

	template <class VD, class ED> CVertex<VD, ED>::CVertex(CGraph<VD, ED>& iGraph, const VD& data): graph(iGraph) {
		index = graph.registerVertex(*this);

		this->data = new VD(data);
	}

	template <class VD, class ED>
	void CVertex<VD, ED>::registerIncomingEdge(CEdge<VD, ED>* e) {
		incomingEdges.push_back(e);
	}

	template <class VD, class ED>
	void CVertex<VD, ED>::registerOutgoingEdge(CEdge<VD, ED>* e) {
		outgoingEdges.push_back(e);
	}

	template <class VD, class ED>
	CEdge<VD, ED>* CVertex<VD, ED>::getIncomingEdge(const CVertex<VD, ED>& predecessorVertex) {
		for (unsigned int q = 0; q < incomingEdges.size(); q++) {
			CEdge<VD, ED>& e = *incomingEdges[q];
			if (&e.getFromVertex() == &predecessorVertex) return &e;
		}

		return NULL;
	}


	template <class VD, class ED>
	vector<CVertex<VD, ED>* > CVertex<VD, ED>::getRelatedVertices(UGraph_Direction graph_direction) {
		CFastSubSet<vector<CVertex<VD, ED>* >, CVertex<VD, ED> > relatedVertices
		= CFastSubSet<vector<CVertex<VD, ED>* >, CVertex<VD, ED> >(graph.getVertices(), true);

		if ((graph_direction == preceding) || (graph_direction == duplex)) {
			const vector<CEdge<VD, ED>* >& incomingEdges = getIncomingEdges();
			for (unsigned int q = 0; q < incomingEdges.size(); q++) {
				relatedVertices.includeElement(incomingEdges[q]->getFromVertex());
			}
		}
		if ((graph_direction == succeeding) || (graph_direction == duplex)) {
			const vector<CEdge<VD, ED>* >& outgoingEdges = getOutgoingEdges();
			for (unsigned int q = 0; q < outgoingEdges.size(); q++) {
				relatedVertices.includeElement(outgoingEdges[q]->getToVertex());
			}
		}

		return relatedVertices.includedsAsVector();
	}

/*	template <class VD, class ED> vector<CVertex<VD, ED>*>* CVertex<VD, ED>::getRelatedVertices(UGraph_Direction graph_direction) {
		int vC = graph.getVertexCount();
		bool vertexIncluded[vC];
		for (int q = 0; q < vC; q++) vertexIncluded[q] = false;

		vector<CVertex<VD, ED>*>* vv = new vector<CVertex<VD, ED>*>();
		vector<CEdge<VD, ED>*>& relatedEdges = getRelatedEdges(graph_direction);
		for (int q = (relatedEdges.size() - 1); q >= 0; q--) {
			CVertex<VD, ED>& fromVertex = relatedEdges[q]->getFromVertex();
			int fromVertexIndex = fromVertex.getIndex();

			if (!vertexIncluded[fromVertexIndex]) {
				vv->push_back(&fromVertex);
				vertexIncluded[fromVertexIndex] = true;
			}
		}

		return vv;
	}*/

	/*template <class VD, class ED> CVertex<VD, ED>& CVertex<VD, ED>::operator = (const CVertex<VD, ED>& other) {
		if (this != &other) {
			graph = other.graph;
			data = other.data;
		}

		return *this;
	};*/
}

#endif /* VERTEX_CPP_ */
