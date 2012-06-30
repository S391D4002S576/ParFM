#ifdef CSETRELATION_H_

#include "Matrix.h"
#include "Module.h"
#include "Cone.h"
#include "Lattice.h"
#include "SetRelation.h"
#include "Flat.h"

#include "Hadron.h"
#include "VectorDescriptor.h"

namespace AlgoTrans {
	template <class R>
	CSetRelation<CFlat<CLattice<R> > > transitionOperation(CSetRelation<CFlat<CLattice<R> > >& a, CSetRelation<CFlat<CLattice<R> > >& b) {
		int aQDim = a.getQDimension();
		int aPDim = a.getPDimension();
		int bQDim = b.getQDimension();
		int bPDim = b.getPDimension();
		int aConstDim = a.getAffineness();
		int bConstDim = b.getAffineness();

		const CMatrix<R>& A = a.matrix();
		const CMatrix<R>& B = b.matrix();

		ASSERT(a.getQDimension() == b.getPDimension());
		int comDim = a.getQDimension();
		ASSERT(a.getAffineness() == b.getAffineness());
		int resultConstDim = a.getAffineness();

		CMatrix<R> (*U)(int dim) = &(CMatrix<R>::getUnitMatrix);
		CMatrix<R> (*Z)(int rows, int cols) = &(CMatrix<R>::getZeroMatrix);

		int aRows = A.getRowCount();
		int bRows = B.getRowCount();

		// Setup Diophantine equations for obtaining the solution
		CMatrix<R> dioph
		= (
			U(aPDim) << Z(aPDim, bQDim + resultConstDim)
				<< -A.getSubMatrix(0, aRows, 0, aPDim).getTransposedMatrix()
				<< Z(aPDim, bRows)
		  ) >> (
			Z(comDim, aPDim + bQDim + resultConstDim)
				<< A.getSubMatrix(0, aRows, aPDim, aQDim).getTransposedMatrix()
				<< -B.getSubMatrix(0, bRows, 0, bPDim).getTransposedMatrix()
		  ) >> (
			Z(bQDim, aPDim) << U(bQDim) << Z(bQDim, resultConstDim)
				<< Z(bQDim, aRows)
				<< -B.getSubMatrix(0, bRows, bPDim, bQDim).getTransposedMatrix()
		  ) >> (
			Z(resultConstDim, aPDim + bQDim) << U(resultConstDim)
				<< -A.getSubMatrix(0, aRows, aPDim + aQDim, aConstDim).getTransposedMatrix()
				<< Z(resultConstDim, bRows)
		  ) >> (
			Z(resultConstDim, aPDim + bQDim) << U(resultConstDim)
				<< Z(resultConstDim, aRows)
				<< -B.getSubMatrix(0, bRows, bPDim + bQDim, bConstDim).getTransposedMatrix()
		  );

		CMatrix<R> diophOrt = dioph.getOrthogonalMatrix();
		CMatrix<R> resultM = diophOrt.getSubMatrixColumns(0, aPDim + bQDim + resultConstDim);
		resultM.reduceToHermiteNormalForm();

		return CSetRelation<CFlat<CLattice<R> > >(aPDim, bQDim,
					CFlat<CLattice<R> >(resultConstDim,
							CLattice<R>(G, matrixToDescriptorSet(resultM))));
	}

	template <class R>
	CSetRelation<CFlat<CModule<R> > > transitionOperation(CSetRelation<CFlat<CModule<R> > >& a, CSetRelation<CFlat<CModule<R> > >& b) {
		int aQDim = a.getQDimension();
		int aPDim = a.getPDimension();
		int bQDim = b.getQDimension();
		int bPDim = b.getPDimension();
		int aAff = a.getAffineness();
		int bAff = b.getAffineness();

		ASSERT(a.getQDimension() == b.getPDimension());
		ASSERT(a.getAffineness() == b.getAffineness()); // This may have to be changed if the relations only carry coefficients for the relevant parameters

		CMatrix<R> A = !(a.matrix());
		CMatrix<R> B = !(b.matrix());
		//A.print();
		//B.print();

		CMatrix<R> (*U)(int dim) = &(CMatrix<R>::getUnitMatrix);
		CMatrix<R> (*Z)(int rows, int cols) = &(CMatrix<R>::getZeroMatrix);

		// Create new modules with common and co-dimensions (a-P ++ common (a-Q / b-P) ++ a-c ++ b-Q ++ b-c)
		CModule<R> aExt = CModule<R>(G, matrixToDescriptorSet(
				(A.getSubMatrixColumns(aAff, aPDim) << -A.getSubMatrixColumns(aAff + aPDim, aQDim)
		         	 << A.getSubMatrixColumns(0, aAff) << Z(A.getRowCount(), bQDim + bAff))
		        >> (Z(bQDim + bAff, aPDim + aQDim + aAff) << U(bQDim + bAff))));

		CModule<R> bExt = CModule<R>(G, matrixToDescriptorSet(
				(Z(B.getRowCount(), aPDim) << B.getSubMatrixColumns(bAff, bPDim)
		                  << Z(B.getRowCount(), aAff) << B.getSubMatrixColumns(bAff + bPDim, bQDim) << B.getSubMatrixColumns(0, bAff))
		        >> (U(aPDim) << Z(aPDim, aQDim + aAff + bQDim + bAff))
			    >> (Z(aAff, aPDim + aQDim) << U(aAff) << Z(aAff, bQDim + bAff))));

		//aExt.print();
		//bExt.print();

		CMatrix<R> intersection = ((aExt * bExt).matrix());

		CMatrix<R> result = (intersection.getSubMatrixColumns(aPDim + aQDim, aAff)
		                        + intersection.getSubMatrixColumns(aPDim + aQDim + aAff + bQDim, bAff))
							<< intersection.getSubMatrixColumns(0, aPDim)
		                    << intersection.getSubMatrixColumns(aPDim + aQDim + aAff, bQDim);

		return CSetRelation<CFlat<CModule<R> > >(aPDim, bQDim,
				CFlat<CModule<R> >(bAff,
				 CModule<R>(C, matrixToDescriptorSet(result))));
	}

	/*template <class R>
	CSetRelation<CFlat<CCone<R> > > transitionOperation(const CSetRelation<CFlat<CCone<R> > >& a, const CSetRelation<CFlat<CCone<R> > >& b) {
		int aQDim = a.getQDimension();
		int aPDim = a.getPDimension();
		int bQDim = b.getQDimension();
		int bPDim = b.getPDimension();
		int aAff = a.getAffineness();
		int bAff = b.getAffineness();

		ASSERT(a.getQDimension() == b.getPDimension());
		ASSERT(a.getAffineness() == b.getAffineness()); // This may have to be changed if the relations only carry coefficients for the relevant parameters

		CMatrix<R> A = a.getConstraintMatrix();
		CMatrix<R> B = b.getConstraintMatrix();

		CMatrix<R> (*U)(int dim) = &(CMatrix<R>::getUnitMatrix);
		CMatrix<R> (*Z)(int rows, int cols) = &(CMatrix<R>::getZeroMatrix);

		// Create new cones with common and co-dimensions (a-P ++ common (a-Q / b-P) ++ a-c ++ b-Q ++ b-c)
		CCone<R> aExt = CCone<R>::fromConstraints(
				(A.getSubMatrixColumns(aAff, aPDim) << -A.getSubMatrixColumns(aAff + aPDim, aQDim)
		         	      << A.getSubMatrixColumns(0, aAff) << Z(A.getRowCount(), bQDim + bAff)));

		CCone<R> bExt = CCone<R>::fromConstraints(
				(Z(B.getRowCount(), aPDim) << B.getSubMatrixColumns(bAff, bPDim)
		                  << Z(B.getRowCount(), aAff) << B.getSubMatrixColumns(bAff + bPDim, bQDim)
		                  << B.getSubMatrixColumns(0, bAff)));

		CCone<R> convexHull = aExt || bExt;

		vector<int> projDims = vector<int>();

		CMatrix<R> result = (intersection.getSubMatrixColumns(aPDim + aQDim, aAff)
		                        + intersection.getSubMatrixColumns(aPDim + aQDim + aAff + bQDim, bAff))
							<< intersection.getSubMatrixColumns(0, aPDim)
		                    << intersection.getSubMatrixColumns(aPDim + aQDim + aAff, bQDim);

		return CSetRelation<CFlat<CCone<R> > >::fromConstraints(aPDim, bQDim, bAff, result);
	}*/
}

#endif
