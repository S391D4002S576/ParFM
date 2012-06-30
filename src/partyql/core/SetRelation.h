#ifndef SETRELATION_H_
#define SETRELATION_H_

#include <string>

#include "Hadron.h"
#include "Descriptor.h"
#include "Domain.h"

namespace AlgoTrans {
	//template <class Y> class SetRelation;
	//template <class Y> SetRelation<Y> operator + (CSetRelation<Y> a, CSetRelation<Y> b);
	//template <class Y> SetRelation<Y> operator * (CSetRelation<Y> a, CSetRelation<Y> b);

	// A relation -- between a pair of sets P and Q of type Y -- specified by a subset of P x Q of type Y
	// ALT: May want to prefer to encapsulate Y instead of inheriting from it (or maybe not)
	template <class Y> class SetRelation: public Y {
		template <class Z> friend bool isRelationSignatureCompatible(SetRelation<Z>& a, SetRelation<Z>& b);
		template <class Z> friend bool operator == (SetRelation<Z> a, SetRelation<Z> b);
		//template <class Z> friend SetRelation<Z> meet(const SetRelation<Z>& a, const SetRelation<Z>& b); // Meet
		//template <class Z> friend SetRelation<Z> join(const SetRelation<Z>& a, const SetRelation<Z>& b); // Join
	private:
		int constDim, pDim, qDim;
	protected:
	public:
		typedef Y DataType;

		void copyRelationSignature(const SetRelation<Y>& other);
		void copyCompatibleRelationSignature(SetRelation<Y>& a, SetRelation<Y>& b);

		//SetRelation() { };
		SetRelation(int iConstDim, int iPDim, int iQDim, const Y& iSet) : Y(iSet), constDim(iConstDim), pDim(iPDim), qDim(iQDim) { };
		SetRelation(const SetRelation<Y>& other) { *this = other; }

		virtual ~SetRelation() { };

		int getConstDimension() const { return constDim; };
		int getPDimension() const { return pDim; };
		int getQDimension() const { return qDim; };

		SetRelation& operator = (const SetRelation& other);

		//SetRelation getSymmetricCounterPart();
		// For relation on one space
		template <class CombinationOperation> SetRelation getReflexiveClosure(const CombinationOperation& combOp);
		SetRelation getTransitiveClosure();
		template <class TransitionCombinationOperation, class CombinationOperation>
		SetRelation getTransitiveClosure_Pugh(const TransitionCombinationOperation& transCombOp, const CombinationOperation& combOp);
		template <class TransitionCombinationOperation, class CombinationOperation>
		SetRelation getReflexiveTransitiveClosure(const TransitionCombinationOperation& transCombOp, const CombinationOperation& combOp) {
			return getReflexiveClosure(combOp).getTransitiveClosure(transCombOp, combOp);
		};

		class Conjunctive_CombinationOperation {
			public:
				static SetRelation evaluate(const SetRelation& a, const SetRelation& b) { return a * b; }
				class UnderlyingOperation : public Y::MeetOperation { };

		};
		class DisjunctiveHull_CombinationOperation {
			public:
				static SetRelation evaluate(const SetRelation& a, const SetRelation& b) { return a + b; }
				class UnderlyingOperation : public Y::MeetOperation { };

		};

		template <class CombinationOperation> class Operations {
		public:
			static SetRelation calcTransition(SetRelation& a, SetRelation& b) { CombinationOperation combOp; return transitionOperation(combOp, a, b); }
			static SetRelation calcCombination(SetRelation& a, SetRelation& b) { return CombinationOperation::evaluate(a, b); }
		};

		typedef Operations<DisjunctiveHull_CombinationOperation> DisjunctiveHullOperations;

		// Text
		std::string toString(bool forceBothDescriptions = false);
		std::string toStringHtml(bool forceBothDescriptions = false);
		std::string toStringJavaScript(bool forceBothDescriptions = false);
		void print(bool forceBothDescriptions = false) { printf("%s\n", toString(forceBothDescriptions).c_str()); }

		static std::string getTypeName() { return (SHORT_TYPENAMES ? "SR<" : "SetRelation<" ) + Y::getTypeName() + ">"; }
		static std::string getTypeNameHtml() { return Html::wrapType((SHORT_TYPENAMES ? "SR" : "SetRelation" ), Y::getTypeNameHtml()); }
	};

	template <class Y> SetRelation<Y> operator + (SetRelation<Y> a, SetRelation<Y> b) { ASSERT(isRelationSignatureCompatible(a, b));
		return SetRelation<Y>(a.getConstDimension(), a.getPDimension(), a.getQDimension(), ((Y) a) + ((Y) b));
	}

	template <class Y> SetRelation<Y> operator || (SetRelation<Y> a, SetRelation<Y> b) { ASSERT(isRelationSignatureCompatible(a, b));
		return SetRelation<Y>(a.getConstDimension(), a.getPDimension(), a.getQDimension(), ((Y) a) || ((Y) b));
	}

	template <class Y> SetRelation<Y> operator * (SetRelation<Y> a, SetRelation<Y> b) {	ASSERT(isRelationSignatureCompatible(a, b));
		return SetRelation<Y>(a.getConstDimension(), a.getPDimension(), a.getQDimension(), ((Y) a) * ((Y) b));
	}

	template <class Y> bool operator == (SetRelation<Y> a, SetRelation<Y> b) {
		return isRelationSignatureCompatible(a, b) && ((Y) a == (Y) b);
	}

	template <class Y> bool operator != (SetRelation<Y> a, SetRelation<Y> b) { return !(a == b); }

	template <class Y> SetRelation<Y>& SetRelation<Y>::operator = (const SetRelation& other) {
		if (this != &other) {
			constDim = other.constDim;
			pDim = other.pDim;
			qDim = other.qDim;
			Y::operator=(other);
		}

		return *this;
	}

	template <class Y>
	bool isRelationSignatureCompatible(SetRelation<Y>& a, SetRelation<Y>& b) {
		return (a.getConstDimension() == b.getConstDimension())
		       && (a.getPDimension() == b.getPDimension())
			   && (a.getQDimension() == b.getQDimension())
			   && (a.getSpaceDimension() == b.getSpaceDimension());
	}

	template <class Y>
	void SetRelation<Y>::copyCompatibleRelationSignature(SetRelation<Y>& a, SetRelation<Y>& b) {
		ASSERT(isRelationSignatureCompatible(a, b) && isSignatureCompatible(a, b));

		copyRelationSignature(a);
	}

	template <class Y>
	void SetRelation<Y>::copyRelationSignature(const SetRelation<Y>& other) {
		setConstDimension(other.getConstDimension());
		setPDimension(other.getPDimension());
		setQDimension(other.getQDimension());

		Y::copySignature(other);
	}

	// a : relation on (P x R), b : relation on (R x Q) --> result: Lamb (p, q) : P x Q . Exi r : R . (p a r) and (r b q)
	template <class CombinationOperation, class Y>
	SetRelation<Y> transitionOperation(const CombinationOperation combOp, SetRelation<Y>& a, SetRelation<Y>& b) {
		ASSERT(a.getQDimension() == b.getPDimension());
		ASSERT(a.getConstDimension() == b.getConstDimension());

		Y aX = ((Y) a ^ Y::universe(b.getQDimension())); // (aC, aP, aQ, [bQ])
		//a.print(); aX.print();
		Y bX = (Y::universe(a.getPDimension()) ^ (Y) b); // ([aP], bC, bP, bQ)
		//b.print(); bX.print();

		// Now permute dimensions so that constants overlap, and result has layout (homog ++ constants ++ p-dim result ++ q-dim result)
		// --> Transform ([aP], bC, bP, bQ) to (bC, [aP], bP, bQ)
		std::vector<int> sourceDims;
		sourceDims.push_back(0);
		for (int q = 0; q < a.getConstDimension(); q++) sourceDims.push_back(1 + a.getPDimension() + q);
		for (int q = 0; q < a.getPDimension(); q++) sourceDims.push_back(1 + q);
		for (int q = 0; q < b.getPDimension() + b.getQDimension(); q++) sourceDims.push_back(1 + a.getPDimension() + a.getConstDimension() + q);
		bX = bX.getDimensionPermutation(sourceDims);
		//bX.print();

		Y abX = CombinationOperation::UnderlyingOperation::evaluate(aX, bX);
		//abX.print();

		// Project out intermediate (aQ / bP) dims
		std::vector<int> retainedDims = std::vector<int>();
		for (int q = 0; q <= a.getConstDimension() + a.getPDimension(); q++) retainedDims.push_back(q);
		for (int q = 0; q < b.getQDimension(); q++) retainedDims.push_back(1 + a.getConstDimension() + a.getPDimension() + a.getQDimension() + q);
		Y abProj = abX.getProjection(retainedDims);

		return SetRelation<Y>(a.getConstDimension(), a.getPDimension(), b.getQDimension(), abProj);
	}

	template <class Y> template <class TransitionCombinationOperation, class CombinationOperation>
	SetRelation<Y> SetRelation<Y>::getTransitiveClosure_Pugh(const TransitionCombinationOperation& transCombOp, const CombinationOperation& combOp) {
		SetRelation<Y> imdR = *this;
		SetRelation<Y> lastR = *this;
		bool first = true;
		while (first || (imdR != lastR)) { first = false;
			imdR.print();
			lastR = imdR;
			SetRelation<Y> transOpResult = transitionOperation(transCombOp, imdR, *this);
			transOpResult.print();
			imdR = combOp.evaluate(imdR, transOpResult);
		}

		return imdR;
	}

	/*template <class Y>
	SetRelation<Y> SetRelation<Y>::getSymmetricCounterPart() {
		Y& org = ((Y) this);
		Description d = org.isDescriptionAvailable(C) ? C : G;

		std::vector<int> sourceDims;
		for (int q = 0; q < getConstDimension() + 1; q++) sourceDims.push(q);
		for (int q = 0; q < getQDimension(); q++) sourceDims.push(1 + getConstDimension() + q + getPDimension());
		for (int q = 0; q < getPDimension(); q++) sourceDims.push(1 + getConstDimension() + q );

		typedef typename Y::DeskriptorSetType DeskriptorSet;
		DeskriptorSet ds = org.getDescriptors(d);

		DeskriptorSet f = ds.getDimensionPermutation(sourceDims);
		Y result =
		return Y();

	}*/

	template <class Y>
	SetRelation<Y> SetRelation<Y>::getTransitiveClosure() {
		ASSERT(getPDimension() == getQDimension());
		typedef typename Y::DeskriptorSetType DeskriptorSet;
		typedef typename DeskriptorSet::DeskriptorType Deskriptor;

		Y a = ((Y) *this).getDual(); // (C, P, Q)
		//a.print(true);
		//a.print(); aX.print();
		//s.print();
		DeskriptorSet eq = DeskriptorSet(1 + getConstDimension() + getPDimension() + getQDimension());
		int sD = getPDimension();
		for (int q = 0; q < sD; q++) {
			eq.addDeskriptor(Deskriptor::getZeroDeskriptor(1 + getConstDimension(), true) << Deskriptor::getUnitDeskriptor(sD, q, true) << Deskriptor::getUnitDeskriptor(sD, q, true));
		}
		for (int q = 0; q < 1 + getConstDimension(); q++) {
			eq.addDeskriptor(Deskriptor::getUnitDeskriptor(1 + getConstDimension() + 2*sD, q, true));
		}
		//Y(C, eq).print(true);
		a = a * Y(C, eq);
		//a.print(true);
		a = a.getDual();
		//a.print(true);
		SetRelation<Y> result = SetRelation<Y>(getConstDimension(), getPDimension(), getQDimension(), a);

		//ASSERT(result = getTransitiveClosure_Pugh()); // Slow test

		return result;
	}
	/* Rethink:
	template <class R>
	CLattice<R> calculateDistanceLattice(SetRelation<CFlat<CLattice<R> > >& rel) {
		int pDim = rel.getPDimension();
		int qDim = rel.getQDimension();
		int cDim = rel.getAffineness();
		const CMatrix<R>& m = rel.matrix();
		ASSERT(pDim == qDim);

		// Short aliases for id- and zero-matrix
		CMatrix<R> (*U)(int dim) = &(CMatrix<R>::getUnitMatrix);
		//CMatrix<R> (*Z)(int rows, int cols) = &(CMatrix<R>::getZeroMatrix);

		CMatrix<R> mTransp = m.getTransposedMatrix();
		CMatrix<R> mTP = mTransp.getSubMatrixRows(cDim, pDim);
		CMatrix<R> mTQ = mTransp.getSubMatrixRows(cDim + pDim, pDim);
		CMatrix<R> mTC = mTransp.getSubMatrixRows(0, cDim);
		CMatrix<R> dioph = U(cDim + pDim) << (mTC >> (mTQ - mTP));

		CLattice<R> result = CLattice<R>(G, matrixToDescriptorSet((!dioph).getSubMatrixColumns(0, pDim + cDim).getHermiteNormalForm()));
		CLattice<R> result2
		= CLattice<R>(G, matrixToDescriptorSet(m.getSubMatrixColumns(0, cDim)
				                      << (m.getSubMatrixColumns(cDim + pDim, pDim) - m.getSubMatrixColumns(cDim, pDim))));

		ASSERT(result == result2);
		//CLattice<R> r3 = rel.getReflexiveTransitiveClosure<CLattice<R> >();
		//ASSERT(result == r3);

		return result;
	}

	Rethink:
	template <class Y> template <class Z>
	Z SetRelation<Y>::getReflexiveTransitiveClosure() {
		ASSERT(getPDimension() == getQDimension());

		return Y::getReflexiveTransitiveClosure();
	}
	*/
}

#include "SetRelation.cpp"
//#include "Transitions.cpp"

#endif /* SETRELATION_H_ */
