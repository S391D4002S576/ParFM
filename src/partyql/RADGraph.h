#ifndef RADGRAPH_H_
#define RADGRAPH_H_

#include <vector>
using std::vector;

#include <list>
using std::list;

#include <algorithm>
using std::min;

#include "graph/Graph.h"
#include "basicmath/Polyheder.h"
#include "basicmath/SetRelation.h"
#include "basicmath/scalar/Interval.h"

namespace AlgoTrans {

	template <class R>
	class CEDGVertexData {
		int statementIndex;
		vector<R> iterationVector;
	};

	template <class R, class Relation>
	class CRADGraph : public CGraph<CPolyheder<R>, Relation > {
		/*CGraph<CEDGVertexData, bool > generateEDG() {
			CGraph<CEDGVertexData, bool > result;

			//vector<vector<int> > dims = vector<vector<int> >();
			vector<vector<CInterval<R> > > vxDimExtremes = vector<vector<CInterval<R> > >();
			for (int v = 0; v < getVertexCount(); v++) {
				CPolyheder<R> poly = getVertex(v).getData();
				vector<CInterval<R> > dimExtremes = vector<CInterval<R> >();
				for (int s = 0; s < poly.getSpaceDimension(); s++) {
					dimExtremes.push_back(poly.asIntervalAlongDimension(s).asIntegralInterval());
				}
				vxDimExtremes.push_back(dimExtremes);
			}

			for (int v = 0; v < getVertexCount(); v++) {
				for (int w = 0; w < getVertexCount(); w++) {

				}
			}
		}*/
	};

	template <class R>
	struct CModuleRADGraph {
		typedef CRADGraph<R, typename CModuleRelation<R>::T > T;
	};

	template <class R>
	struct CLatticeRADGraph {
		typedef CRADGraph<R, typename CLatticeRelation<R>::T > T;
	};


	template <class R, class Relation>
	struct CRADGVertex {
		typedef CVertex<CPolyheder<R>, Relation > T;
	};

	template <class R>
	struct CModuleRADGVertex {
		typedef typename CRADGVertex<R, typename CModuleRelation<R>::T >::T T;
	};

	template <class R>
	struct CLatticeRADGVertex {
		typedef typename CRADGVertex<R, typename CLatticeRelation<R>::T >::T T;
	};


	template <class R, class Relation>
	struct CRADGEdge {
		typedef CEdge<CPolyheder<R>, Relation > T;
	};

	template <class R>
	struct CModuleRADGEdge {
		typedef typename CRADGEdge<R, typename CModuleRelation<R>::T >::T T;
	};

	template <class R>
	struct CLatticeRADGEdge {
		typedef typename CRADGEdge<R, typename CLatticeRelation<R>::T >::T T;
	};
}

#endif /* RADGRAPH_H_ */
