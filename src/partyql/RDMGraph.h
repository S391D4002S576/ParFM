/* OBSOLETE */

#ifndef RDMGRAPH_H_
#define RDMGRAPH_H_

#include <vector>
using std::vector;

#include <list>
using std::list;

#include <algorithm>
using std::min;

#include <sstream>
using std::ostream;

#include "graph/Graph.h"
#include "basicmath/AffineTransformation.h"
#include "basicmath/Polyheder.h"
#include "basicmath/Matrix.h" 
#include "RADGraph.h"
#include "basicmath/SetRelation.h"

#include "SpacePartitioner.h"

#define PAPER1
#ifdef PAPER1
#include "../utils/SVGELattice.h"
#include "../utils/SVGElement.h"
#endif

namespace AlgoTrans {
	template <class R> class CRDMGraph: public CGraph<CPolyheder<R>, CAffineTransformation<R> > {
		typedef CVertex<CPolyheder<R>, CAffineTransformation<R> > CRDMGVertex;
		typedef CEdge<CPolyheder<R>, CAffineTransformation<R> > CRDMGEdge;
	private:		
	protected:
		vector<CRDMGVertex* > computationVertices;
		vector<CRDMGVertex* > dataVertices;
		 
	public:
		#ifdef PAPER1
			// Verbose output for paper
//			ostream* osLatex;
	//		ostream* osHtml;
			string verboseOutputPath;
			string verboseOutputPrefix;
			CVector<int> latticeLowerBounds;
			CVector<int> latticeUpperBounds;
			int latticeBasicWidth;
			vector<CSVGElement* > resultingSVGLattices;
			vector<CLattice<R>* > resultingLattices;
		#endif
		
		
		int parmCount; // XXX: We should use a more advanced mechanism for this..
		
		CRDMGraph() {  };
		
		CRDMGraph(const CRDMGraph<R>& iOriginal) { *this = iOriginal; };
		
		CRDMGraph& operator = (const CRDMGraph& other) {
			if (this != &other) {
				parmCount = other.parmCount;
				computationVertices = other.computationVertices;
				dataVertices = other.dataVertices;
				CGraph<CPolyheder<R>, CAffineTransformation<R> >::operator=(other);
			}
			
			return *this;
		}

		int getComputationVertexCount() { return computationVertices.size(); }
		int getDataVertexCount() { return dataVertices.size(); }
		
		CRDMGVertex& getComputationVertex(int index) { return *computationVertices[index]; }
		CRDMGVertex& getDataVertex(int index) { return *dataVertices[index]; }
		
		CRDMGVertex& addComputationVertex(const CPolyheder<R>& data) {
			CRDMGVertex& v = addVertex(data);
			
			computationVertices.push_back(&v);
			
			return v;
		}

		CRDMGVertex& addDataVertex(const CPolyheder<R>& data) {
			CRDMGVertex& v = addVertex(data);
			
			dataVertices.push_back(&v);
			
			return v;
		}
		
		typename CModuleRADGraph<R>::T calculateAffinelyApproximatedRSDG();		
		typename CLatticeRADGraph<R>::T calculateFlatLatticeApproximatedRSDG();
		
		CompGraphSpacePartitioner<R>* performLatticeBasedSpacePartitioning();
		CompGraphSpacePartitioner<R>* performFlatBasedSpacePartitioning();
		
		string toStringHtml();
	};
}

#include "RDMGraph.cpp"

#endif /* RDMGRAPH_H_ */
