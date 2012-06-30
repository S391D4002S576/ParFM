#ifdef PARTITIONRELATION_H_

#include "scalar/Integer.h"

namespace AlgoTrans {
	template <class Y> string CPartitionRelation<Y>::toString() const {
		string result = "PartitionRelation(" + CInteger(dims.size()).toString() + ", ";

		for (unsigned int q = 0; q < dims.size(); q++) {
			result += CInteger(dims[q]).toString();
			result += ", ";
		}

		result += Y::toString() + ")";

		return result;
	}

	template <class Y> template <class R>
	CPartitionRelation<Y> CPartitionRelation<Y>::getPairRelationFromCoupledConstraints(int leftDim, int rightDim, const CMatrix<R>& leftConstraints, const CMatrix<R>& rightConstraints) {
		CPartitionRelation<Y> result = CPartitionRelation<Y>(leftDim + rightDim);

		result.setByConstraints(leftConstraints << rightConstraints);

		return result;
	}

	template <class Y> template <class R>
	CPartitionRelation<Y> CPartitionRelation<Y>::getPairRelationFromIntersectedConstraints(int leftDim, int rightDim, const CMatrix<R>& leftConstraints, const CMatrix<R>& rightConstraints) {
		CPartitionRelation<Y> result = CPartitionRelation<Y>(leftDim + rightDim);

		CMatrix<R> (*Z)(int rows, int cols) = &(CMatrix<R>::getZeroMatrix);

		result.setByConstraints((leftConstraints << Z(leftConstraints.getRowCount(), rightDim))
								>> (Z(rightConstraints.getRowCount(), leftDim) <<  rightConstraints));

		return result;
	}

	template <class R>
	CPolyheder<R> getStrictlyPositiveFunctionPolyheder(int dim, CPolyheder<R> polyheder) {
		CVector<R> (*UV)(int dim, int unitDim) = &(CVector<R>::getUnitVector);
		CVector<R> (*ZV)(int dim) = &(CVector<R>::getZeroVector);

		CMatrix<R> consM = polyheder.getFarkasPolyheder().getConstraintMatrix();

		//consM.print();
		CMatrix<R> resultFuncs = CMatrix<R>(1 + dim);
		for (int q = 0; q < consM.getRowCount(); q++) {
			resultFuncs.addRow(ZV(1) << consM[q]);
			resultFuncs.addRow(UV(1, 0) << (consM[q] + UV(dim, 0)));
		}

		return CPolyheder<R>::fromInequalities(CCone<R>::fromGenerators(resultFuncs).getConstraintMatrix());
	}

	template <class R>
	CPolyheder<R> getReflexiveValidPartitions(int affineness, int spaceDim, CPolyheder<R> polyheder) {
		CMatrix<R> (*Z)(int rows, int cols) = &(CMatrix<R>::getZeroMatrix);
		CMatrix<R> (*U)(int dim) = &(CMatrix<R>::getUnitMatrix);
		CVector<R> (*UV)(int dim, int unitDim) = &(CVector<R>::getUnitVector);
		CVector<R> (*ZV)(int dim) = &(CVector<R>::getZeroVector);

		CCone<R> dualGen = CCone<R>::fromGenerators(polyheder.getFarkasPolyheder().getConstraintMatrix());
		CMatrix<R> constZ = (U(affineness) << Z(affineness, 2*spaceDim));
		CMatrix<R> coeffEq = Z(spaceDim, affineness) << U(spaceDim) << U(spaceDim);
		CMatrix<R> posDualGen = (dualGen && CCone<R>::fromConstraints(coeffEq >> -coeffEq
	              								>> constZ >> -constZ)).getGeneratorMatrix();
		CMatrix<R> funs = posDualGen.getSubMatrixColumns(0, affineness)
							<< posDualGen.getSubMatrixColumns(affineness + spaceDim, spaceDim);

		CMatrix<R> funDis = CMatrix<R>(1 + affineness + spaceDim);
		for (int q = 0; q < funs.getRowCount(); q++) {
			funDis.addRow(ZV(2) << funs[q].getSubVector(1, affineness - 1 + spaceDim));
			funDis.addRow(UV(2, 0) << funs[q].getSubVector(1, affineness - 1 + spaceDim).getVectorDividedByGCDOfElements());
		}
		funDis = funDis >> (Z(affineness, 1) << U(affineness) << Z(affineness, spaceDim))
						>> -(Z(affineness, 1) << U(affineness) << Z(affineness, spaceDim));

		return CPolyheder<R>::fromInequalities(CCone<R>::fromGenerators(funDis).getConstraintMatrix());
	}

	template <class R>
	CPolyheder<R> getReflexiveValid2Partitions(int affineness, int spaceDim, CPolyheder<R> polyheder) {
		CMatrix<R> (*Z)(int rows, int cols) = &(CMatrix<R>::getZeroMatrix);
		CMatrix<R> (*U)(int dim) = &(CMatrix<R>::getUnitMatrix);
		//CVector<R> (*UV)(int dim, int unitDim) = &(CVector<R>::getUnitVector);
		//CVector<R> (*ZV)(int dim) = &(CVector<R>::getZeroVector);

		CMatrix<R> consM = polyheder.getFarkasPolyheder().getConstraintMatrix();
		CMatrix<R> superConeM
		= (consM << Z(consM.getRowCount(), spaceDim))
		  >>  (Z(spaceDim, affineness) << -U(spaceDim) << U(spaceDim) << -U(spaceDim))
		  >> -(Z(spaceDim, affineness) << -U(spaceDim) << U(spaceDim) << -U(spaceDim));
		CCone<R> superCone = CCone<R>(C, matrixToDescriptorTeam(superConeM));

		for (int q = 0; q < 2*spaceDim; q++) {
			superCone = superCone.getProjectionFourierMotzkin(affineness, q == 0);
		}
		for (int q = 0; q < affineness; q++) {
			superCone = superCone.getProjectionFourierMotzkin(0, false);
		}
		CMatrix<R> funs = superCone.matrix();
		CMatrix<R> funDis = Z(funs.getRowCount(), affineness) << funs;

		return CPolyheder<R>::fromInequalities(funDis);
	}

	template <class R>
	CPolyheder<R> getReflexiveDismissingPartitions(int affineness, int spaceDim, CPolyheder<R> polyheder) {
		CMatrix<R> (*Z)(int rows, int cols) = &(CMatrix<R>::getZeroMatrix);
		CMatrix<R> (*U)(int dim) = &(CMatrix<R>::getUnitMatrix);
		CVector<R> (*UV)(int dim, int unitDim) = &(CVector<R>::getUnitVector);
		CVector<R> (*ZV)(int dim) = &(CVector<R>::getZeroVector);

		//polyheder.getFarkasPolyheder().getConstraintMatrix().print();
		CCone<R> dualGen = CCone<R>::fromGenerators(polyheder.getFarkasPolyheder().getConstraintMatrix());
		//dualGen.print();
		//(!dualGen).print();
		CMatrix<R> constZ = (U(affineness) << Z(affineness, 2*spaceDim)).getSubMatrixRows(1, affineness - 1);
		CMatrix<R> coeffEq = Z(spaceDim, affineness) << U(spaceDim) << U(spaceDim);
		CMatrix<R> posDualGen = (dualGen && CCone<R>::fromConstraints(coeffEq >> -coeffEq
	              								>> -UV(2*spaceDim + affineness, 0) // const must be negative for dismissing function
	              								>> constZ >> -constZ)).getGeneratorMatrix();
		//posDualGen.print();
		CMatrix<R> funs = posDualGen.getSubMatrixColumns(0, affineness)
							<< posDualGen.getSubMatrixColumns(affineness + spaceDim, spaceDim);
		//funs.print();

		CMatrix<R> funDis = CMatrix<R>(1 +  affineness + spaceDim);
		for (int q = 0; q < funs.getRowCount(); q++) {
			funDis.addRow(ZV(2) << funs[q].getSubVector(1, affineness - 1 + spaceDim));
			if (!funs[q][0].isZero()) {
				funDis.addRow(UV(2, 0) << funs[q].getSubVector(1, affineness - 1 + spaceDim).getVectorDividedByGCDOfElements());
			}
		}
		funDis = funDis >> (Z(affineness, 1) << U(affineness) << Z(affineness, spaceDim))
						>> -(Z(affineness, 1) << U(affineness) << Z(affineness, spaceDim));
		//funDis.print();

		return CPolyheder<R>::fromInequalities(CCone<R>::fromGenerators(funDis).getConstraintMatrix());
	}

	template <class Y> template <class R>
	CPartitionRelation<Y>
	CPartitionRelation<Y>::getDismissingRelationFromPolyhedralRelation(vector<int> dims, int affineness, int leftIx, int rightIx, CPolyheder<R> polyheder) {
		CMatrix<R> (*Z)(int rows, int cols) = &(CMatrix<R>::getZeroMatrix);
		CMatrix<R> (*U)(int dim) = &(CMatrix<R>::getUnitMatrix);
		//CVector<R> (*UV)(int dim, int unitDim) = &(CVector<R>::getUnitVector);

		int lowIx = (leftIx < rightIx) ? leftIx : rightIx;
		int upIx = (leftIx > rightIx) ? leftIx : rightIx;

		int z = 0;
		for (int q = 0; q < lowIx; q++) z += dims[q];
		int leftBorderDim = z;
		z = 0;
		for (int q = lowIx + 1; q < upIx; q++) z += dims[q];
		int innerBandDim = z;
		z = 0;
		for (unsigned int q = upIx + 1; q < dims.size(); q++) z += dims[q];
		int rightBorderDim = z;

		if (leftIx == rightIx) {
			int d = dims[lowIx] - affineness;
			CPolyheder<R> dismissingFuncs = getReflexiveDismissingPartitions(affineness, d, polyheder);
			//dismissingFuncs.print(); CCone<R>::fromConstraints(dismissingFuncs.getFarkasPolyheder().getConstraintMatrix()).getGeneratorMatrix().print();

			ASSERT(!dismissingFuncs.isEmpty());
			CMatrix<R> genM = dismissingFuncs.getFarkasPolyheder().getConstraintMatrix();
			CMatrix<R> resultM = Z(genM.getRowCount(), leftBorderDim) << genM << Z(genM.getRowCount(), rightBorderDim);

			return CPartitionRelation<Y>::fromGenerators(dims, resultM);
		} else {
			CMatrix<R> constrM = CCone<R>::fromConstraints(getStrictlyPositiveFunctionPolyheder(dims[leftIx] + dims[rightIx] - affineness, polyheder).getFarkasPolyheder().getConstraintMatrix()).getGeneratorMatrix();

			CMatrix<R> leftM = constrM.getSubMatrixColumns(0, affineness) << constrM.getSubMatrixColumns(affineness, dims[leftIx] - affineness);
			CMatrix<R> rightM = Z(constrM.getRowCount(), affineness) << -constrM.getSubMatrixColumns(dims[leftIx], dims[rightIx] - affineness);

			int rowc = constrM.getRowCount();
			CMatrix<R> resultM =
			  (Z(rowc, leftBorderDim)
			  << ((lowIx == leftIx) ? leftM : rightM)   << Z(rowc, innerBandDim)
			  << ((lowIx == leftIx) ? rightM : leftM)   << Z(rowc, rightBorderDim))
			  >>
			  (Z(affineness, leftBorderDim)
			  << U(affineness) << Z(affineness, (dims[lowIx] - affineness) + innerBandDim)
			  << -U(affineness) << Z(affineness, (dims[upIx] - affineness) + rightBorderDim))
			  ;

			return CPartitionRelation<Y>::fromConstraints(dims, resultM);
		}
	}

	template <class Y> template <class R>
	CPartitionRelation<Y>
	CPartitionRelation<Y>::getValidRelationFromPolyhedralRelation(vector<int> dims, int affineness, int leftIx, int rightIx, CPolyheder<R> polyheder) {
		CMatrix<R> (*Z)(int rows, int cols) = &(CMatrix<R>::getZeroMatrix);
		CMatrix<R> (*U)(int dim) = &(CMatrix<R>::getUnitMatrix);
		//CVector<R> (*UV)(int dim, int unitDim) = &(CVector<R>::getUnitVector);

		CMatrix<R> constrM = polyheder.getConstraintMatrix();

		int lowIx = (leftIx < rightIx) ? leftIx : rightIx;
		int upIx = (leftIx > rightIx) ? leftIx : rightIx;

		int z = 0;
		for (int q = 0; q < lowIx; q++) z += dims[q];
		int leftBorderDim = z;
		z = 0;
		for (int q = lowIx + 1; q < upIx; q++) z += dims[q];
		int innerBandDim = z;
		z = 0;
		for (unsigned int q = upIx + 1; q < dims.size(); q++) z += dims[q];
		int rightBorderDim = z;

		if (leftIx == rightIx) {
			int d = dims[lowIx] - affineness;
			//CPolyheder<R> validFuncs = getReflexiveValidPartitions(affineness, d, polyheder);
			//polyheder.print();
			//validFuncs.print();
			CPolyheder<R> validFuncs = getReflexiveValid2Partitions(affineness, d, polyheder);
			//validFuncs2.print();
			//ASSERT(validFuncs == validFuncs2);
			//validFuncs.print(); CCone<R>::fromConstraints(validFuncs.getFarkasPolyheder().getConstraintMatrix()).getGeneratorMatrix().print();

			//ASSERT(!validFuncs.isEmpty());
			CMatrix<R> genM = validFuncs.getFarkasPolyheder().getConstraintMatrix();
			CMatrix<R> resultM = Z(genM.getRowCount(), leftBorderDim) << genM << Z(genM.getRowCount(), rightBorderDim);

			return CPartitionRelation<Y>(dims, Y(C, matrixToDescriptorTeam(resultM)));
		} else {
			CMatrix<R> leftM = constrM.getSubMatrixColumns(0, affineness) << constrM.getSubMatrixColumns(affineness, dims[leftIx] - affineness);
			CMatrix<R> rightM = Z(constrM.getRowCount(), affineness) << -constrM.getSubMatrixColumns(dims[leftIx], dims[rightIx] - affineness);

			int rowc = constrM.getRowCount();
			CMatrix<R> resultM =
			  (Z(rowc, leftBorderDim)
			  << ((lowIx == leftIx) ? leftM : rightM)   << Z(rowc, innerBandDim)
			  << ((lowIx == leftIx) ? rightM : leftM)   << Z(rowc, rightBorderDim))
			  >>
			  (Z(affineness, leftBorderDim)
			  << U(affineness) << Z(affineness, (dims[lowIx] - affineness) + innerBandDim)
			  << -U(affineness) << Z(affineness, (dims[upIx] - affineness) + rightBorderDim))
			  ;

			return CPartitionRelation<Y>(dims, Y(C, matrixToDescriptorTeam(resultM)));
		}
	}


/*	template <class Y> template <class R>
	CPartitionRelation<Y>
	CPartitionRelation<Y>::getValid2RelationFromPolyhedralRelation(vector<int> dims, int affineness, int leftIx, int rightIx, CPolyheder<R> polyheder) {
		CMatrix<R> (*Z)(int rows, int cols) = &(CMatrix<R>::getZeroMatrix);
		CMatrix<R> (*U)(int dim) = &(CMatrix<R>::getUnitMatrix);
		//CVector<R> (*UV)(int dim, int unitDim) = &(CVector<R>::getUnitVector);

		CMatrix<R> constrM = polyheder.getConstraintMatrix();

		int lowIx = (leftIx < rightIx) ? leftIx : rightIx;
		int upIx = (leftIx > rightIx) ? leftIx : rightIx;

		int z = 0;
		for (int q = 0; q < lowIx; q++) z += dims[q];
		int leftBorderDim = z;
		z = 0;
		for (int q = lowIx + 1; q < upIx; q++) z += dims[q];
		int innerBandDim = z;
		z = 0;
		for (unsigned int q = upIx + 1; q < dims.size(); q++) z += dims[q];
		int rightBorderDim = z;

		if (leftIx == rightIx) {
			int d = dims[lowIx] - affineness;
			CPolyheder<R> validFuncs = getReflexiveValid2Partitions(affineness, d, polyheder);
			//validFuncs.print(); CCone<R>::fromConstraints(validFuncs.getFarkasPolyheder().getConstraintMatrix()).getGeneratorMatrix().print();

			ASSERT(!validFuncs.isEmpty());
			CMatrix<R> genM = validFuncs.getFarkasPolyheder().getConstraintMatrix();
			CMatrix<R> resultM = Z(genM.getRowCount(), leftBorderDim) << genM << Z(genM.getRowCount(), rightBorderDim);

			return CPartitionRelation<Y>::fromGenerators(dims, resultM);
		} else {
			CMatrix<R> leftM = constrM.getSubMatrixColumns(0, affineness) << constrM.getSubMatrixColumns(affineness, dims[leftIx] - affineness);
			CMatrix<R> rightM = Z(constrM.getRowCount(), affineness) << -constrM.getSubMatrixColumns(dims[leftIx], dims[rightIx] - affineness);

			int rowc = constrM.getRowCount();
			CMatrix<R> resultM =
			  (Z(rowc, leftBorderDim)
			  << ((lowIx == leftIx) ? leftM : rightM)   << Z(rowc, innerBandDim)
			  << ((lowIx == leftIx) ? rightM : leftM)   << Z(rowc, rightBorderDim))
			  >>
			  (Z(affineness, leftBorderDim)
			  << U(affineness) << Z(affineness, (dims[lowIx] - affineness) + innerBandDim)
			  << -U(affineness) << Z(affineness, (dims[upIx] - affineness) + rightBorderDim))
			  ;

			return CPartitionRelation<Y>::fromConstraints(dims, resultM);
		}
	}*/

	template <class Y>
	CPartitionRelation<Y> meet(const CPartitionRelation<Y>& a, const CPartitionRelation<Y>& b) {
		CPartitionRelation<Y> result = CPartitionRelation<Y>();
		result.copyCompatibleRelationSignature(a, b);

		result.setByConstraints(a.getConstraintMatrix() >> b.getConstraintMatrix());

		return result;
	}

	template <class Y> CPartitionRelation<Y> operator &&(const CPartitionRelation<Y>& a, const CPartitionRelation<Y>& b) {
		return meet(a, b);
	}

	template <class Y>
	CPartitionRelation<Y> CPartitionRelation<Y>::source(int iSpaceDimension) {
		return CPartitionRelation(iSpaceDimension, Y::source(iSpaceDimension));
	}

	template <class Y>
	CPartitionRelation<Y> CPartitionRelation<Y>::universe(int iSpaceDimension) {
		return CPartitionRelation(iSpaceDimension, Y::universe(iSpaceDimension));
	}

	template <class Y>
	void CPartitionRelation<Y>::copyCompatibleRelationSignature(const CPartitionRelation<Y>& a, const CPartitionRelation<Y>& b) {
		dims.clear();
		vector<int> aDims = a.getDimGroupDims();
		vector<int> bDims = b.getDimGroupDims();
		ASSERT(aDims.size() == bDims.size());
		for (unsigned int q = 0; q < aDims.size(); q++) {
			dims.push_back(aDims[q]);
			ASSERT(aDims[q] == bDims[q]);
		}
	}
}

#endif
