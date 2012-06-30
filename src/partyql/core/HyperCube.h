#ifndef LATTICE_H_
#define LATTICE_H_

#include <vector>
using std::vector;

#include <list>
using std::list;

#include <algorithm>
using std::min;

#include "Vektor.h"
#include "Descriptor.h"
#include "Hadron.h"

namespace AlgoTrans {
	template <class R> class CHyperCube : public CHadron<CBidirectionalVectorDescriptorSet<R> > {
		typedef CBidirectionalVectorDescriptorSet<R> DescriptorSet;

		//template <class G> friend CHyperCube<G> fromGenerators(const CMatrix<G>& generatorMatrix);

		//template <class G> friend class CHermiteNormalFormReducer;
		template <class G> friend bool operator == (CHyperCube<G> a, CHyperCube<G> b);
		//template <class G> friend CHyperCube<G> intersect(const CHyperCube<G>& a, const CHyperCube<G>& b);
		template <class G> friend CHyperCube<G> operator * (CHyperCube<G> a, CHyperCube<G> b); // intersection of generators
		template <class G> friend CHyperCube<G> operator + (CHyperCube<G> a, CHyperCube<G> b); // union of generators
		template <class G> friend CHyperCube<G> operator ! (CHyperCube<G> a); // Orthogonal HyperCube
		//template <class G> friend void calculateMinkowskiSum(CHyperCube<G>& result, const CHyperCube<G>& a, const CHyperCube<G>& b);
		//template <class G> friend void calculateIntersection(CHyperCube<G>& result, const CHyperCube<G>& a, const CHyperCube<G>& b);

	private:
		//CHyperCube* getProjectedHyperCube(CVector<int> retainedDims, CVector<R> projectDimValues);
	protected:
		CHyperCube() { };
	public:
		CHyperCube<R> getProjection(vector<int> retainedDimensiones) const;
		//CHyperCube<R> getReflexiveTransitiveClosure(int spaceDim);

		CHyperCube(int iSpaceDimension, HadronDescription description) : CHadron<DescriptorSet>(iSpaceDimension, description) { }
		CHyperCube(const CHyperCube<R>& iOriginal) : CHadron<DescriptorSet>(iOriginal) { }
		CHyperCube(HadronDescription description, const DescriptorSet& descriptors) : CHadron<DescriptorSet>(description, descriptors) { };
		//CHyperCube(const CHadron<DescriptorSet>& iOriginal) : CHadron<DescriptorSet>(iOriginal) { };

		void orthogonalize();
		CHyperCube<R> getOrthogonalHyperCube() const;

		//string toStringLatex() const;

		//CHyperCube<R>& operator = (const CHyperCube<R>& other);


		// Affine stuff
		// Find an exemplar vector of the affine HyperCube represented by this homogenised HyperCube
		//CVector<R>* getAffineExemplar() const;
		//bool containsAffinePoint(const CVector<R>& point);
		//bool containsPoint(const CVector<R>& point);

		//bool isEmptyAffine() const;

		//CVector<R> getAffineExemplar() const;
		//CHyperCube<R> getLinearHyperCube() const;
	};

	template <class R> CHyperCube<R> operator + (CHyperCube<R> a, CHyperCube<R> b) {
		return CHyperCube<R>(((CHadron<CBidirectionalVectorDescriptorSet<R> >) a) + ((CHadron<CBidirectionalVectorDescriptorSet<R> >) b));
	}

	template <class R> CHyperCube<R> operator * (CHyperCube<R> a, CHyperCube<R> b) {
		return CHyperCube<R>(((CHadron<CBidirectionalVectorDescriptorSet<R> >) a) * ((CHadron<CBidirectionalVectorDescriptorSet<R> >) b));
	}

	template <class R> CHyperCube<R> operator !(CHyperCube<R>& a) {
		ASSERTM(false, "Not yet implemented");
		return CHyperCube<R>(!((CHadron<CBidirectionalVectorDescriptorSet<R> >) a));
	}

	template <class R> bool operator == (CHyperCube<R> a, CHyperCube<R> b) {
		// Convert both sets of generating vertices to HNF and check whether they are equal
		ASSERTM(false, "Not yet implemented");
		// Minimize twice + sort + compare element-wise
		// Or simplify a xor b ?
		/*CMatrix<R> xA = a.matrix().getHermiteNormalForm();
		CMatrix<R> xB = b.matrix().getHermiteNormalForm();

		int commDim = min(xA.getRowCount(), xB.getRowCount());

		for (int q = 0; q < commDim; q++) if (!(xA[q] == xB[q])) return false;
		for (int q = commDim; q < xA.getRowCount(); q++) if (!(xA[q].isZero())) return false;
		for (int q = commDim; q < xB.getRowCount(); q++) if (!(xB[q].isZero())) return false;*/

		return true;
	}


	/*template <class R> bool CHyperCube<R>::containsPoint(const CVector<R>& point) {
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
	}*/


	/*template <class R> CHyperCube<R> CHyperCube<R>::getReflexiveTransitiveClosure(int spaceDim) {
		int aff = this->getDimension() - 2*spaceDim;

		CMatrix<R> m = matrix();
		CMatrix<R> result
		= m.getSubMatrixColumns(1, aff - 1)
		  << (m.getSubMatrixColumns(aff + spaceDim, spaceDim) - m.getSubMatrixColumns(aff, spaceDim));

		return CHyperCube<R>(G, matrixToDescriptorSet(result.getHermiteNormalForm()));
	}

	template <class R>
	const CMatrix<R> CHyperCube<R>::matrix() {
		return DescriptorSetToMatrix(CHadron<CBidirectionalVectorDescriptorSet<R> >::getDescriptors(G));
	}*/
}

#endif /* HyperCube_H_ */
