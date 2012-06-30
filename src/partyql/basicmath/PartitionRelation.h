#ifndef PARTITIONRELATION_H_
#define PARTITIONRELATION_H_

#include <vector>
using std::vector;

#include <list>
using std::list;

#include <algorithm>
using std::min;

#include "Matrix.h"
#include "Hadron.h"

namespace AlgoTrans {
	template <class Y> 
	class CPartitionRelation: public Y {
		//template <class Z> friend bool isRelationSignatureCompatible(const CSetRelation<Z>& a, const CSetRelation<Z>& b);
		//template <class Z> friend bool operator == (const CSetRelation<Z>& a, const CSetRelation<Z>& b);
		template <class Z> friend CPartitionRelation<Z> meet(const CPartitionRelation<Z>& a, const CPartitionRelation<Z>& b); // Meet
		template <class Z> friend CPartitionRelation<Z> join(const CPartitionRelation<Z>& a, const CPartitionRelation<Z>& b); // Join

	private:
		vector<int> dims;
 
	public:
		CPartitionRelation(int iGlobalDim, HadronDescription description) : Y(iGlobalDim, description) { };
		CPartitionRelation(vector<int> iDims, const Y& iY) : Y(iY), dims(iDims) { }
		/*CPartitionRelation(vector<int> iDims) {
			for (unsigned int q = 0; q < iDims.size(); q++) dims.push_back(iDims[q]);
		};*/
		
		/*template <class R> static CPartitionRelation<Y>
		fromConstraints(vector<int> dims, const CMatrix<R>& iConstraintMatrix);
		template <class R> static CPartitionRelation<Y>
		fromGenerators(vector<int> dims, const CMatrix<R>& iGeneratorMatrix);
		
		template <class R> static CPartitionRelation<Y> source(int iSpaceDimension);
		template <class R> static CPartitionRelation<Y> universe(int iSpaceDimension);*/
		
		static CPartitionRelation<Y> source(int iSpaceDimension);
		static CPartitionRelation<Y> universe(int iSpaceDimension);
				
		void copyCompatibleRelationSignature(const CPartitionRelation<Y>& a, const CPartitionRelation<Y>& b);
		
		void addDimGroup(int dimension) { dims.push_back(dimension); }
		vector<int> getDimGroupDims() const { return dims; }
		
		template <class R>
		static CPartitionRelation<Y> 
		getPairRelationFromCoupledConstraints(int leftDim, int rightDim, const CMatrix<R>& leftConstraints, const CMatrix<R>& rightConstraints); 
				
		template <class R>
		static CPartitionRelation<Y> 
		getPairRelationFromIntersectedConstraints(int leftDim, int rightDim, const CMatrix<R>& leftConstraints, const CMatrix<R>& rightConstraints);
		
		template <class R>
		static CPartitionRelation<Y>
		getValidRelationFromPolyhedralRelation(vector<int> dims, int affineness, int leftIx, int rightIx, CPolyheder<R> polyheder);
		template <class R>
		static CPartitionRelation<Y>
		getDismissingRelationFromPolyhedralRelation(vector<int> dims, int affineness, int leftIx, int rightIx, CPolyheder<R> polyheder);
		//CCone<R> getProjectionFourierMotzkin(int projectionDim) const;
		
		//CCone<R> getSimplifiedConeFast() const;
		
		string toString() const;		
//		string toStringHtml() const { return ("Cone(" + constraintMatrix.toStringHtml() + ")"); };
		void print() const { printf((toString() + "\n").c_str()); };
		
		//CCone<R>& operator = (const CCone<R>& other);
	};

	template <class Y> CPartitionRelation<Y> meet(const CPartitionRelation<Y>& a, const CPartitionRelation<Y>& b);
	template <class Y> CPartitionRelation<Y> join(const CPartitionRelation<Y>& a, const CPartitionRelation<Y>& b);
	
	template <class Y> 
	CPartitionRelation<Y> operator * (CPartitionRelation<Y> a, CPartitionRelation<Y> b) {
		return CPartitionRelation<Y>((Y) a * (Y) b);
	}
	
	template <class Y> 
	CPartitionRelation<Y> operator + (CPartitionRelation<Y> a, CPartitionRelation<Y> b) {
		return CPartitionRelation<Y>((Y) a + (Y) b);
	}

}

#include "PartitionRelation.cpp"

#endif /*CONE_H_*/
