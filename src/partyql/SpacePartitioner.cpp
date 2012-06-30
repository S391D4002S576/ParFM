/* OBSOLETE */

#ifdef SPACEPARTITIONER_H_

namespace AlgoTrans {
	template <class R> CompGraphSpacePartitioner<R>::CompGraphSpacePartitioner(CRDMGraph<R>* iRDMGraph) {
		rdmGraph = iRDMGraph;
	}

	template <class R> void CompGraphSpacePartitioner<R>::performAffineLatticeSpacePartitioning() {
		// Convert polyhedra representation to reduced graph of partitioning relations
		if (osHtml != NULL) *osHtml << "<b>Input (Polyhedral representation)</b><br/>" << rdmGraph->toStringHtml() << "<hr/>";
		
		rdmGraph->osHtml = osHtml;
		RADLatticeGraph flApproxRSDG = rdmGraph->calculateFlatLatticeApproximatedRSDG();
			
		flApproxRSDG.osHtml = osHtml;
		dependenceLatticeClosure = (RADLatticeGraph*) flApproxRSDG.calcJoinTransitiveClosure();

		// Determine distance lattices 
		if (osHtml != NULL) *osHtml << "<hr/><b>" << "Calculation of distance lattices:" << "</b>" << "<br/>";
		for (int q = 0; q < dependenceLatticeClosure->getVertexCount(); q++) {
			CRADGLatticeVertex& v = (*dependenceLatticeClosure)(q);
			LatticeRelation l = LatticeRelation(dependenceLatticeClosure->getUniqueEdgeBetween(v, v).getData());
			distanceLattices.push_back(new CLattice<R>(calculateDistanceLattice(l)));
			if (osHtml != NULL ) *osHtml << "Statement " << q << ":<br/>"
			<< distanceLattices[q]->toStringHtml() << "<br/>";
		}
		
		// Polyhedral decomposition
		decomposedGraph = new CRDMGraph<R>();
		for (int p = 0; p < rdmGraph->getComputationVertexCount(); p++) { CRDMGVertex& vxP = rdmGraph->getComputationVertex(p);
			CRDMGVertex& newVXP = decomposedGraph->addComputationVertex(vxP.getData().polyhedralLatticeDecomposition(*distanceLattices[p]));
			for (int d = 0; d < rdmGraph->getDataVertexCount(); d++) { CRDMGVertex& vxD = rdmGraph->getDataVertex(d);
				if (p == 0) decomposedGraph->addDataVertex(vxD.getData());
				vector<CRDMGEdge* > edges = rdmGraph->getEdgesBetween(vxP, vxD);
				for (unsigned int e = 0; e < edges.size(); e++) { CRDMGEdge& egE = *edges[e];
					CAffineTransformation<R> at = egE.getData().polyhedralLatticeDecomposition(*distanceLattices[p]);
					decomposedGraph->addEdge(newVXP, decomposedGraph->getDataVertex(d), at);
				}
			}				
		}
		
		if (osHtml != NULL) *osHtml << "<hr/><b>Output (Polyhedral representation)</b><br/>" << decomposedGraph->toStringHtml(); 
	}
	
	template <class R> void CompGraphSpacePartitioner<R>::performAffineSetSpacePartitioning() {
		// Convert polyhedra representation to reduced graph of partitioning relations
		if (osHtml != NULL) *osHtml << "<b>Input (Polyhedral representation)</b><br/>" << rdmGraph->toStringHtml() << "<hr/>";
		
		rdmGraph->osHtml = osHtml;
		RADModuleGraph flApproxRSDG = rdmGraph->calculateAffinelyApproximatedRSDG();
			
		flApproxRSDG.osHtml = osHtml; 
		dependenceFlatClosure = (RADModuleGraph*) flApproxRSDG.calcJoinTransitiveClosure();
	}
}

#endif /* SPACEPARTITIONER_H_ */
