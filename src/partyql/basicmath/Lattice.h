#ifndef LATTICE_H_
#define LATTICE_H_

#include <vector>
using std::vector;

#include <list>
using std::list;

#include <algorithm>
using std::min;

#include "Vektor.h"
#include "Matrix.h"
#include "factorization/HermiteNormalFormReduction.h"
#include "Hadron.h"

#include "VectorDescriptor.h"

namespace AlgoTrans {
	template <class R> class CLattice : public CHadron<CBidirectionalVectorDescriptorSet<R> > {
		typedef CBidirectionalVectorDescriptorSet<R> DescriptorSet;

		template <class G> friend CLattice<G> fromGenerators(const CMatrix<G>& generatorMatrix);

		template <class G> friend class CHermiteNormalFormReducer;
		template <class G> friend bool operator == (CLattice<G> a, CLattice<G> b);
		template <class G> friend CLattice<G> intersect(const CLattice<G>& a, const CLattice<G>& b);
		template <class G> friend CLattice<G> operator * (CLattice<G> a, CLattice<G> b); // Intersection
		template <class G> friend CLattice<G> operator + (CLattice<G> a, CLattice<G> b); // Minkowski sum
		template <class G> friend CLattice<G> operator ! (const CLattice<G>& a); // Orthogonal Lattice
		template <class G> friend void calculateMinkowskiSum(CLattice<G>& result, const CLattice<G>& a, const CLattice<G>& b);
		template <class G> friend void calculateIntersection(CLattice<G>& result, const CLattice<G>& a, const CLattice<G>& b);

	private:
		CLattice* getProjectedLattice(CVector<int> retainedDims, CVector<R> projectDimValues);
	protected:
		CLattice() { };
	public:
		const CMatrix<R> matrix();
		//CLattice& operator >>= (const CLattice<R>& m);
		//CLattice& operator >>= (const CVector<R>& v);

		vector<R> getPolyhedralDecompositionMultipliers();

		CLattice<R> getProjection(vector<int> retainedDimensiones) const;
		CLattice<R> getReflexiveTransitiveClosure(int spaceDim);

		CLattice(int iSpaceDimension, HadronDescription description) : CHadron<DescriptorSet>(iSpaceDimension, description) { }
		CLattice(const CLattice<R>& iOriginal) : CHadron<DescriptorSet>(iOriginal) { }
		CLattice(HadronDescription description, const DescriptorSet& descriptors) : CHadron<DescriptorSet>(description, descriptors) { };
		CLattice(const CHadron<DescriptorSet>& iOriginal) : CHadron<DescriptorSet>(iOriginal) { };

		//void addGeneratingVertex(CVector<R>* v);
		//void addGeneratingVertex(const CVector<R>& v);


		void orthogonalize();
		CLattice<R> getOrthogonalLattice() const;

		//int getGeneratingVectorCount() const { return matrix.getRowCount(); };

		//string toStringLatex() const;

		//CLattice<R>& operator = (const CLattice<R>& other);

		//int getAffineSpaceDimension() const { return getDimension() - 1; };

		bool isEmpty() const;

		// Affine stuff
		// Find an exemplar vector of the affine lattice represented by this homogenised lattice
		//CVector<R>* getAffineExemplar() const;
		bool containsAffinePoint(const CVector<R>& point);
		bool containsPoint(const CVector<R>& point);

		bool isEmptyAffine() const;

		CVector<R> getAffineExemplar() const;
		CLattice<R> getLinearLattice() const;
	};

	template <class R> CLattice<R> operator + (CLattice<R> a, CLattice<R> b) {
		return CLattice<R>(((CHadron<CBidirectionalVectorDescriptorSet<R> >) a) + ((CHadron<CBidirectionalVectorDescriptorSet<R> >) b));
	}

	template <class R> CLattice<R> operator * (CLattice<R> a, CLattice<R> b) {
		return CLattice<R>(((CHadron<CBidirectionalVectorDescriptorSet<R> >) a) * ((CHadron<CBidirectionalVectorDescriptorSet<R> >) b));
	}

	template <class R> CLattice<R> meet(CLattice<R>& a, CLattice<R>& b) {
		return CLattice<R>(meet((CHadron<CBidirectionalVectorDescriptorSet<R> >&) a), ((CHadron<CBidirectionalVectorDescriptorSet<R> >&) b));
	}

	template <class R> CLattice<R> join(CLattice<R>& a, CLattice<R>& b) {
		return CLattice<R>(join((CHadron<CBidirectionalVectorDescriptorSet<R> >&) a), ((CHadron<CBidirectionalVectorDescriptorSet<R> >&) b));
	}

	template <class R> CLattice<R> operator !(CLattice<R>& a) {
		return CLattice<R>(!((CHadron<CBidirectionalVectorDescriptorSet<R> >&) a));
	}

	template <class R> vector<R> CLattice<R>::getPolyhedralDecompositionMultipliers() {
		vector<R> multipliers;

		CMatrix<R> coreMatrix = CMatrix<R>(matrix);
		CHermiteNormalFormReducer<R> reducer = CHermiteNormalFormReducer<R>(coreMatrix);

		reducer.reduce();

		for (int p = 0; p < reducer.pivotColumns.size(); p++) {
			multipliers.push_back(coreMatrix[p][reducer.pivotColumns]);
		}

		return multipliers;
	}

	template <class R> bool operator == (CLattice<R> a, CLattice<R> b) {
		// Convert both sets of generating vertices to HNF and check whether they are equal
		CMatrix<R> xA = a.matrix().getHermiteNormalForm();
		CMatrix<R> xB = b.matrix().getHermiteNormalForm();

		int commDim = min(xA.getRowCount(), xB.getRowCount());

		for (int q = 0; q < commDim; q++) if (!(xA[q] == xB[q])) return false;
		for (int q = commDim; q < xA.getRowCount(); q++) if (!(xA[q].isZero())) return false;
		for (int q = commDim; q < xB.getRowCount(); q++) if (!(xB[q].isZero())) return false;

		return true;
	}

	/*template <class R> void CLattice<R>::addGeneratingVertex(CVector<R>* v) {
		matrix.addRow(v);
	}

	template <class R> void CLattice<R>::addGeneratingVertex(const CVector<R>& v) {
		matrix.addRow(v);
	}

	template <class R> CLattice<R>* CLattice<R>::getUnitMatrixLattice(int dim) {
		return new CLattice<R>(CMatrix<R>::getUnitMatrix(dim));
	}*/

	/*template <class R> bool CLattice<R>::containsAffinePoint(const CVector<R>& vector) const {
		CMatrix<R> homogV = vector << U(1);

		CMatrix<R> thisHNF = getHermiteNormalForm();
		CMatrix<R> thisHNFHomogV = (*this >> homogV).getHermiteNormalForm();

		return (thisHNF.getGeneratingVectorCount() == thisHNFHomogV.getGeneratingVectorCount());
	}*/

	template <class R> CVector<R> CLattice<R>::getAffineExemplar() const {
		CMatrix<R> m
		= matrix.getSubMatrixColumns(matrix.getColumnCount() - 1, 1)
		  << matrix.getSubMatrixColumns(0, matrix.getColumnCount() - 1);
		m.reduceToHermiteNormalForm();

		ASSERT(m[0][0] == 1);

		return m[0].getSubVector(1, m.getColumnCount() - 1);
	}

	template <class R> CLattice<R> CLattice<R>::getLinearLattice() const {
		CMatrix<R> m
		= matrix.getSubMatrixColumns(matrix.getColumnCount() - 1, 1)
		  << matrix.getSubMatrixColumns(0, matrix.getColumnCount() - 1);
		m.reduceToHermiteNormalForm();

		ASSERT(m[0][0] == 1);

		return CLattice<R>::fromGenerators(m.getSubMatrix(1, m.getRowCount() - 1, 1, m.getColumnCount() - 1));
	}

	/*template <class R> bool CLattice<R>::isEmpty() const {
		return (getGeneratingVectorCount() == 0);
	}

	template <class R> bool CLattice<R>::isEmptyAffine() const {
		// Checks whether we can find an exemplar vector
		int hCoord = getGeneratingVectorCount() - 1;
		R homog = (*this)[getSpaceDimension() - 1][hCoord];

		// Find a generating vector G for which GCD(homogeneousCoord(G), homog) == 1
		CVector<R> hResult;
		for (int q = 0; q < getGeneratingVectorCount(); q++) {
			CExtendedGCDResults<R> xGCDResults = extendedGCD(homog, (*this)[q][hCoord]);
			if (xGCDResults.gcd.isOne()) return true;
		}

		return false;
	}*/

	template <class R> bool CLattice<R>::containsAffinePoint(const CVector<R>& point) {
		CMatrix<R> mTranspose = DescriptorSetToMatrix(this->getDescriptors(G)).getTransposedMatrix();

		// Short aliases for id- and zero-matrix
		CMatrix<R> (*U)(int dim) = &(CMatrix<R>::getUnitMatrix);
		CMatrix<R> (*Z)(int rows, int cols) = &(CMatrix<R>::getZeroMatrix);

		int aDim = this->getDimension() - 1;

		CMatrix<R> pointMatrixT = ((CMatrix<R>) point).getTransposedMatrix();

		CMatrix<R> diopheq
		=    (Z(aDim, 1) << U(aDim)    << -mTranspose.getSubMatrixRows(1, aDim))
		  >> (-U(1) << Z(1, aDim) << mTranspose.getSubMatrixRows(0, 1))
		  >> (-pointMatrixT << U(aDim)    << Z(aDim, mTranspose.getColumnCount()));

		CMatrix<R> orthM = diopheq.getOrthogonalMatrix().getHermiteNormalForm();
		if (orthM.getRowCount() == 0) return false;

		R hGcd = orthM[0][0];
		for (int q = 1; q < orthM.getRowCount(); q++) {
			if (hGcd == 1) return true;
			hGcd = gcd(hGcd, orthM[q][0]);
		}
		return (hGcd.getAbs() == 1);
	}

	template <class R> bool CLattice<R>::containsPoint(const CVector<R>& point) {
		CMatrix<R> mTranspose = DescriptorSetToMatrix(this->getDescriptors(G)).getTransposedMatrix();

		// Short aliases for id- and zero-matrix
		CMatrix<R> (*U)(int dim) = &(CMatrix<R>::getUnitMatrix);
		CMatrix<R> (*Z)(int rows, int cols) = &(CMatrix<R>::getZeroMatrix);

		int aDim = this->getDimension();

		CMatrix<R> pointMatrix = CMatrix<R>(point.getLength());
		pointMatrix.addRow(point);

		CMatrix<R> diopheq
		=    (U(aDim)    << -mTranspose.getSubMatrixRows(0, aDim) << Z(aDim, 1))
		  >> (U(aDim)    << Z(aDim, mTranspose.getColumnCount()) << -pointMatrix.getTransposedMatrix());

		CMatrix<R> orthM = !diopheq;
		if (orthM.getRowCount() == 0) return false;

		R hGcd = orthM[0][orthM.getColumnCount() - 1];
		for (int q = 1; q < orthM.getRowCount(); q++) {
			if (hGcd == 1) return true;
			hGcd = gcd(hGcd, orthM[q][orthM.getColumnCount() - 1]);
		}
		return (hGcd.getAbs() == 1);
	}

	template <class R> CLattice<R>* CLattice<R>::getProjectedLattice(CVector<int> retainedDims, CVector<R> projDimValues) {
		CMatrix<R> dioph = CMatrix<R>(retainedDims.getLength() + matrix().getRowCount() + 1);

		// Short aliases for id- and zero-matrix
		//CMatrix<R> (*U)(int dim) = &(CMatrix<R>::getUnitMatrix);
		//CMatrix<R> (*Z)(int rows, int cols) = &(CMatrix<R>::getZeroMatrix);
		CVector<R> (*UV)(int dim, int unitDim) = &(CVector<R>::getUnitVector);
		CVector<R> (*ZV)(int dim) = &(CVector<R>::getZeroVector);

		CMatrix<R> thisTransp = matrix().getTransposedMatrix();

		int nextDim = 0;
		int p = 0;
		int r = 0;
		for (int q = 0; q < retainedDims.getLength(); q++) { int d = retainedDims[q];
			while (nextDim < d) {
				dioph >>= (ZV(retainedDims.getLength()) << thisTransp[nextDim++] << -projDimValues[p++]);
			}
			dioph >>= (UV(retainedDims.getLength(), r++) << -thisTransp[nextDim++] << ZV(1));
		}

		while (nextDim < this->getDimension()) {
			dioph >>= (ZV(retainedDims.getLength()) << thisTransp[nextDim++] << -projDimValues[p++]);
		}

		CMatrix<R> sols = (!dioph).getSubMatrixColumns(0, retainedDims.getLength()).getHermiteNormalForm();

		return new CLattice<R>(sols);
	}

	template <class R> CLattice<R> CLattice<R>::getProjection(vector<int> retainedDimensions) const {
		std::sort(retainedDimensions.begin(), retainedDimensions.end());

		CMatrix<R> matrix = matrix();
		CMatrix<R> result;
		for (unsigned int q = 0; q < retainedDimensions.size(); q++) {
			ASSERT((0 <= retainedDimensions[q]) && (retainedDimensions[q] < matrix.getColumnCount()));
			result = result << matrix.getSubMatrixColumns(retainedDimensions[q], 1);
		}

		return CLattice<R>(matrixToDescriptorSet(result.getHermiteNormalForm()), G);
	}

	template <class R> CLattice<R> CLattice<R>::getReflexiveTransitiveClosure(int spaceDim) {
		int aff = this->getDimension() - 2*spaceDim;

		CMatrix<R> m = matrix();
		CMatrix<R> result
		= m.getSubMatrixColumns(1, aff - 1)
		  << (m.getSubMatrixColumns(aff + spaceDim, spaceDim) - m.getSubMatrixColumns(aff, spaceDim));

		return CLattice<R>(G, matrixToDescriptorSet(result.getHermiteNormalForm()));
	}

	template <class R>
	const CMatrix<R> CLattice<R>::matrix() {
		return DescriptorSetToMatrix(CHadron<CBidirectionalVectorDescriptorSet<R> >::getDescriptors(G));
	}
}

#endif /* LATTICE_H_ */
