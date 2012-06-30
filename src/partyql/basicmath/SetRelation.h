#ifndef CSETRELATION_H_
#define CSETRELATION_H_

#include <string>
using std::string;

#include "Matrix.h"
#include "Module.h"
#include "Lattice.h"
#include "Flat.h"

#include "Hadron.h"
#include "Descriptor.h"

namespace AlgoTrans {
	template <class Y> class CSetRelation;
	template <class Y> CSetRelation<Y> operator + (CSetRelation<Y> a, CSetRelation<Y> b);
	template <class Y> CSetRelation<Y> operator * (CSetRelation<Y> a, CSetRelation<Y> b);

	/* A relation
	 * between a pair of sets P and Q of type Y
	 * specified by a subset of P x Q of type Y
	 */
	// XXX: May want to prefer to encapsulate Y instead of inheriting from it (or maybe not)
	template <class Y> class CSetRelation: public Y {
		template <class Z> friend bool isRelationSignatureCompatible(CSetRelation<Z>& a, CSetRelation<Z>& b);
		template <class Z> friend bool operator == (CSetRelation<Z> a, CSetRelation<Z> b);
		template <class Z> friend CSetRelation<Z> meet(const CSetRelation<Z>& a, const CSetRelation<Z>& b); // Meet
		template <class Z> friend CSetRelation<Z> join(const CSetRelation<Z>& a, const CSetRelation<Z>& b); // Join
	private:
		int pDim, qDim;
	protected:
	public:
		// SuperSpace dim is not adjusted by these methods!
		//void setPDimension(int newPDimension) { pDim = newPDimension; }
		//void setQDimension(int newQDimension) { qDim = newQDimension; }

		void copyRelationSignature(const CSetRelation<Y>& other);
		void copyCompatibleRelationSignature(CSetRelation<Y>& a, CSetRelation<Y>& b);

		//CSetRelation() { };
		CSetRelation(int iPDim, int iQDim, const Y& iSet) : Y(iSet), pDim(iPDim), qDim(iQDim) { };
		CSetRelation(const CSetRelation<Y>& other);

		/*template <class R>
		static CSetRelation<Y> fromConstraints(int iPDim, int iQDim, int iConstDim, const CMatrix<R>& constraintMatrix);
		template <class R>
		static CSetRelation<Y> fromGenerators(int iPDim, int iQDim, int iConstDim, const CMatrix<R>& generatorMatrix);
		*/

		virtual ~CSetRelation() { };

		int getPDimension() const { return pDim; }
		int getQDimension() const { return qDim; }

		CSetRelation& operator = (const CSetRelation& other);

		string toString(bool forceBothDescriptions = false);
		string toStringHtml(bool forceBothDescriptions = false);
		void print(bool forceBothDescriptions = false);

		template <class Z>
		Z getReflexiveTransitiveClosure();

		class ConjunctiveOperations {
		public:
			static CSetRelation calcTransition(CSetRelation& a, CSetRelation& b) { return transitionOperation(a, b); }
			static CSetRelation calcCombination(CSetRelation& a, CSetRelation& b) { return a * b; }
		};

		class DisjunctiveHullOperations {
		public:
			static CSetRelation calcTransition(CSetRelation& a, CSetRelation& b) { return transitionOperation(a, b); }
			static CSetRelation calcCombination(CSetRelation& a, CSetRelation& b) { return a + b; }
		};
	};

	template <class R> struct CLatticeRelation {
		typedef CSetRelation<CFlat<CLattice<R> > > T;
	};

	template <class R> struct CModuleRelation {
		typedef CSetRelation<CFlat<CModule<R> > > T;
	};


	template <class Y> CSetRelation<Y> operator + (CSetRelation<Y> a, CSetRelation<Y> b) {
		ASSERT(isRelationSignatureCompatible(a, b));

		return CSetRelation<Y>(a.getPDimension(), a.getQDimension(), ((Y) a) + ((Y) b));
	}

	template <class Y> CSetRelation<Y> operator * (CSetRelation<Y> a, CSetRelation<Y> b) {
		ASSERT(isRelationSignatureCompatible(a, b));

		return CSetRelation<Y>(a.getPDimension(), a.getQDimension(), ((Y) a) * ((Y) b));
	}


}

#include "SetRelation.cpp"
#include "Transitions.cpp"

#endif /* SETRELATION_H_ */
