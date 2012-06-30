/* OBSOLETE */

#ifndef SPACEPARTITIONER_H_
#define SPACEPARTITIONER_H_

#include "basicmath/Lattice.h"

#include "RADGraph.h"

#include "basicmath/Polyheder.h"
#include "basicmath/AffineTransformation.h"
#include "basicmath/SetRelation.h"

#include <vector>
using std::vector;

#include <sstream>
using std::ostream;

namespace AlgoTrans {
	template <class R> class CRDMGraph;
	
	template <class R>
	class CompGraphSpacePartitioner {
		typedef CVertex<CPolyheder<R>, CAffineTransformation<R> > CRDMGVertex;
		typedef CEdge<CPolyheder<R>, CAffineTransformation<R> > CRDMGEdge;
		typedef typename CLatticeRADGraph<R>::T RADLatticeGraph;
		typedef typename CModuleRADGraph<R>::T RADModuleGraph;
		typedef typename CLatticeRelation<R>::T LatticeRelation;
		typedef typename CRADGVertex<R, LatticeRelation>::T CRADGLatticeVertex;
		typedef typename CRADGEdge<R, LatticeRelation>::T CRADGLatticeEdge;
		
		protected:
			CRDMGraph<R>* rdmGraph;
			vector<CLattice<R>* > distanceLattices;
			
			CRDMGraph<R>* decomposedGraph;
			
		public:
			RADLatticeGraph* dependenceLatticeClosure;
			RADModuleGraph* dependenceFlatClosure;
			
			ostream* osHtml;
			
			CompGraphSpacePartitioner(CRDMGraph<R>* iRDMGraph);
			
			CAffineTransformation<R> polyhedrallyDecomposeAffineTransformation(CAffineTransformation<R>& at, CLattice<R>& core) {
				return CAffineTransformation<R>(polyhedrallyDecomposeConstraintMatrix(CMatrix<R>(at.getMatrix()), core));
			}
			
			void performAffineLatticeSpacePartitioning();			
			void performAffineSetSpacePartitioning();
	};
}

#include "SpacePartitioner.cpp"

#endif /*SPACEPARTITIONER_H_*/
