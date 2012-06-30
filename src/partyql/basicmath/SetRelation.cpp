#ifndef CSETRELATION_CPP_
#define CSETRELATION_CPP_

#include "SetRelation.h"
#include "scalar/Integer.h"
#include "Flat.h"

#include "Hadron.h"
#include "VectorDescriptor.h"

#include "../../utils/Html.h"

namespace AlgoTrans {
	template <class Y> CSetRelation<Y>::CSetRelation(const CSetRelation<Y>& other)
	{
		*this = other;
	}

	/*template <class Y>
	CSetRelation<Y>(int iPDim, int iQDim, int iConstDim, Description description, const DescriptorTeam<R>& descriptors)
	: Y(description, descriptors), pDim(iPDim), qDim(iQDim), constDim(iConstDim) {
	}*/

	/*template <class Y> template <class R>
	CSetRelation<Y> CSetRelation<Y>::fromGenerators(int iPDim, int iQDim, int iConstDim, const CMatrix<R>& generatorMatrix) {
		CSetRelation<Y> result = CSetRelation<Y>(iPDim, iQDim, iConstDim);

		result.setByGenerators(generatorMatrix);

		return result;
	}*/

	template <class Y> CSetRelation<Y>& CSetRelation<Y>::operator = (const CSetRelation& other) {
		if (this != &other) {
			pDim = other.pDim;
			qDim = other.qDim;
			Y::operator=(other);
		}

		return *this;
	}

	template <class Y> bool operator == (CSetRelation<Y> a, CSetRelation<Y> b) {
		return isRelationSignatureCompatible(a, b) && ((Y) a == (Y) b);
	}

	template <class Y> bool operator != (CSetRelation<Y> a, CSetRelation<Y> b) {
		return !(a == b);
	}

	template <class Y> string CSetRelation<Y>::toString(bool forceBothDescriptions) {
		return "SetRelation(" + CInteger(getPDimension()).toString() + ", " + CInteger(getQDimension()).toString() + ", "
		       + Y::toString(forceBothDescriptions) + ")";
	}

	template <class Y> string CSetRelation<Y>::toStringHtml(bool forceBothDescriptions) {
		return Html::wrapper("<table><tr><td>SetRelation</td><td>P-dim: " + CInteger(getPDimension()).toString() + "</td><td>Q-dim: " + CInteger(getQDimension()).toString() + "</td></tr></table>",
		       Y::toStringHtml(forceBothDescriptions), std::string("BBCCAA"));
	}

	template <class Y> void CSetRelation<Y>::print(bool forceBothDescriptions) {
		printf("%s\n", toString(forceBothDescriptions).c_str());
	}

	template <class Y>
	bool isRelationSignatureCompatible(CSetRelation<Y>& a, CSetRelation<Y>& b) {
		return (a.getPDimension() == b.getPDimension())
			   && (a.getQDimension() == b.getQDimension())
			   && isSignatureCompatible(a, b);
	}

	template <class Y>
	void CSetRelation<Y>::copyCompatibleRelationSignature(CSetRelation<Y>& a, CSetRelation<Y>& b) {
		ASSERT(isRelationSignatureCompatible(a, b) && isSignatureCompatible(a, b));

		copyRelationSignature(a);
	}

	template <class Y>
	void CSetRelation<Y>::copyRelationSignature(const CSetRelation<Y>& other) {
		setPDimension(other.getPDimension());
		setQDimension(other.getQDimension());

		Y::copySignature(other);
	}

	/*
	template <class Y>
	CSetRelation<Y> meet(const CSetRelation<Y >& a, const CSetRelation<Y >& b) {
		ASSERT(isRelationSignatureCompatible(a, b));

		CSetRelation<Y >& result = CSetRelation<Y >(a.getPDimension(), a.getQDimension());

		Y::result = ((Y) a) && ((Y) b);

		return result;
	}

	template <class R>
	CSetRelation<Y> join(const CSetRelation<Y>& a, const CSetRelation<Y>& b) {
		ASSERT(isRelationSignatureCompatible(a, b));

		CSetRelation<Y >& result = CSetRelation<Y >(a.getPDimension(), a.getQDimension());

		Y::result = ((Y) a) + ((Y) b);

		return result;
	}*/

	template <class R>
	CLattice<R> calculateDistanceLattice(CSetRelation<CFlat<CLattice<R> > >& rel) {
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

	template <class Y> template <class Z>
	Z CSetRelation<Y>::getReflexiveTransitiveClosure() {
		ASSERT(getPDimension() == getQDimension());

		return Y::getReflexiveTransitiveClosure();
	}
}

#endif /* SETRELATION_CPP_ */
