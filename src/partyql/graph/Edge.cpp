#ifndef EDGE_CPP_
#define EDGE_CPP_

#include "Graph.h"

namespace AlgoTrans {
	template <class VD, class ED>  
	CEdge<VD, ED>::CEdge(CGraph<VD, ED>& iGraph, CVertex<VD, ED>& iFromVertex, CVertex<VD, ED>& iToVertex) :
		graph(iGraph), fromVertex(iFromVertex), toVertex(iToVertex) {
		graph.registerEdge(this);
		
		data = NULL;
	}
	
	template <class VD, class ED> 
	CEdge<VD, ED>::CEdge(CGraph<VD, ED>& iGraph, CVertex<VD, ED>& iFromVertex, CVertex<VD, ED>& iToVertex,
								 const ED& iData) :
		graph(iGraph), fromVertex(iFromVertex), toVertex(iToVertex) {
		graph.registerEdge(this);
		
		data = new ED(iData);
	}
	
	template <class VD, class ED> 
	CEdge<VD, ED>::CEdge(CGraph<VD, ED>& iGraph, CVertex<VD, ED>& iFromVertex, CVertex<VD, ED>& iToVertex,
								 ED* iData) :
		graph(iGraph), fromVertex(iFromVertex), toVertex(iToVertex) {
		graph.registerEdge(this);
		
		data = iData;
	}
}

#endif /* EDGE_CPP_ */
