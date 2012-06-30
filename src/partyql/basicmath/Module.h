#ifndef MODULE_H_
#define MODULE_H_

#include <vector>
using std::vector;

#include <list>
using std::list;

#include <algorithm>
using std::min;

#include "Vektor.h"
#include "VectorDescriptor.h"
#include "Matrix.h"
#include "Hadron.h"
#include "factorization/HermiteNormalFormReduction.h"

namespace AlgoTrans {
	template <class R> class CModule : public CHadron<CBidirectionalVectorDescriptorSet<R> > {
		typedef CBidirectionalVectorDescriptorSet<R> DescriptorSet;

		template <class G> friend class CHermiteNormalFormReducer;
		template <class G> friend bool operator == (CModule<G> a, CModule<G> b);
		//template <class G> friend CModule<G> intersect(const CModule<G>& a, const CModule<G>& b);
		//template <class G> friend CModule<G> operator * (CModule<G>& a, CModule<G>& b);
		//template <class G> friend CModule<G> operator + (CModule<G>& a, CModule<G>& b);
		//template <class G> friend CModule<G> operator ! (CModule<G>& a);
		//template <class G> friend void calculateMinkowskiSum(CModule<G>& result, const CModule<G>& a, const CModule<G>& b);
		//template <class G> friend void calculateIntersection(CModule<G>& result, const CModule<G>& a, const CModule<G>& b);

	protected:
		CModule() { };
	public:
		const CMatrix<R> matrix();

		CModule(int iSpaceDimension, HadronDescription description) : CHadron<DescriptorSet>(iSpaceDimension, description) { }
		CModule(const CModule<R>& iOriginal) : CHadron<DescriptorSet>(iOriginal) {  };
		CModule(HadronDescription description, const DescriptorSet& descriptors) : CHadron<DescriptorSet>(description, descriptors) { };
		CModule(const CHadron<DescriptorSet>& iOriginal) : CHadron<DescriptorSet>(iOriginal) { };

		//void setSpaceDimension(int newSpaceDimension);

		//static CModule<R>* getUnitMatrixModule(int dim);
		//void addGeneratingVertex(CVector<R>* v);
		//CModule<R> getOrthogonalModule() const;

		// Projection represented in lower-dimensional (projected) space
		CModule<R> getProjection(vector<int> retainedDimensions); // obsolete v
		CModule<R> getReflexiveTransitiveClosure(int spaceDim);

		//int getSpaceDimension() const { return matrix.getColumnCount(); };
		//int getGeneratingVectorCount() const { return matrix.getRowCount(); };

		//string toString() const;
		//string toStringHtml() const { return matrix.toStringHtml(); };
		//void print() const;

		//CModule<R>& operator = (const CModule<R>& other);

		bool isEmpty();
		bool isAffineEmpty();
	};

	template <class R> CModule<R> operator + (CModule<R> a, CModule<R> b) {
		return CModule<R>(((CHadron<CBidirectionalVectorDescriptorSet<R> >) a) + ((CHadron<CBidirectionalVectorDescriptorSet<R> >) b));
	}

	template <class R> CModule<R> operator * (CModule<R> a, CModule<R> b) {
		return CModule<R>(((CHadron<CBidirectionalVectorDescriptorSet<R> >) a) * ((CHadron<CBidirectionalVectorDescriptorSet<R> >) b));
	}

	template <class R> CModule<R> meet(CModule<R>& a, CModule<R>& b) {
		return CModule<R>(meet((CHadron<CBidirectionalVectorDescriptorSet<R> >&) a), ((CHadron<CBidirectionalVectorDescriptorSet<R> >&) b));
	}

	template <class R> CModule<R> join(CModule<R>& a, CModule<R>& b) {
		return CModule<R>(join((CHadron<CBidirectionalVectorDescriptorSet<R> >&) a), ((CHadron<CBidirectionalVectorDescriptorSet<R> >&) b));
	}

	template <class R> CModule<R> operator !(CModule<R> a) {
		return CModule<R>(!((CHadron<CBidirectionalVectorDescriptorSet<R> >) a));
	}

	template <class R> bool operator == (CModule<R> a, CModule<R> b) {
		// Convert both sets of generating vertices to HNF and check whether they are equal
		CMatrix<R> xA = a.matrix().getHermiteNormalForm();
		CMatrix<R> xB = b.matrix().getHermiteNormalForm();

		int commDim = min(xA.getRowCount(), xB.getRowCount());

		for (int q = 0; q < commDim; q++) if (!(xA[q] == xB[q])) return false;
		for (int q = commDim; q < xA.getRowCount(); q++) if (!(xA[q].isZero())) return false;
		for (int q = commDim; q < xB.getRowCount(); q++) if (!(xB[q].isZero())) return false;

		return true;
	}

	template <class R> CModule<R> CModule<R>::getProjection(vector<int> retainedDimensions) {
		std::sort(retainedDimensions.begin(), retainedDimensions.end());

		CMatrix<R> result = CMatrix<R>(0);
		for (unsigned int q = 0; q < retainedDimensions.size(); q++) {
			ASSERT((0 <= retainedDimensions[q]) && (retainedDimensions[q] < matrix().getColumnCount()));
			result = result << matrix().getSubMatrixColumns(retainedDimensions[q], 1);
		}

		return CModule<R>(G, matrixToDescriptorSet(result.getHermiteNormalForm()));
	}

	template <class R> CModule<R> CModule<R>::getReflexiveTransitiveClosure(int spaceDim) {
		int aff = this->getDimension() - 2*spaceDim;
		CMatrix<R> result
		= matrix().getSubMatrixColumns(1, aff - 1)
		  << (matrix().getSubMatrixColumns(aff + spaceDim, spaceDim) - matrix().getSubMatrixColumns(aff, spaceDim));

		return CModule<R>(G, matrixToDescriptorSet(result.getHermiteNormalForm()));
	}

	template <class R>
	const CMatrix<R> CModule<R>::matrix() {
		return DescriptorSetToMatrix(CHadron<CBidirectionalVectorDescriptorSet<R> >::getDescriptors(G));
	}
}

#endif /*MODULE_H_*/
