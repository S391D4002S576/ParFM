#ifndef HERMITENORMALFORMREDUCTION_CPP_
#define HERMITENORMALFORMREDUCTION_CPP_

#include "HermiteNormalFormReduction.h"

namespace AlgoTrans {
	template <class R> CHermiteNormalFormReducer<R>::CHermiteNormalFormReducer(CMatrix<R>& iMatrix): X(iMatrix) {

	}

	template <class R> void CHermiteNormalFormReducer<R>::reduce() {
		rank = 0;

		for (int j = 0; j < X.getColumnCount(); j++) {
			bool pivotFound = false;
			R pivot;
			for (int i = rank; i < X.getRowCount(); i++) {
				if (X(i, j) != 0) {
					if (!pivotFound) { // Found it now! :)
						X.swapVertices(rank, i);
						pivotFound = true;
						pivot = X(rank, j);
						pivotColumns.push_back(j);

						// We must ensure it's positive
						if (pivot < 0) {
							X[rank] = -X[rank];
							pivot = -pivot;
						}
					} else {
						CExtendedExtendedGCDResults<R> xxGCD = extendedExtendedGCD(pivot, X(i, j));

						X.combineVertices(rank, i, xxGCD.leftBezoutFactor, xxGCD.rightBezoutFactor, xxGCD.leftLCMCofactor, -xxGCD.rightLCMCofactor);

						pivot = xxGCD.gcd;
					}
				}
			}

			// Now reduce the rows above the rank-row
			if (pivotFound) {
				for (int i = rank - 1; i >= 0; i--) {
					R quotient = X(i, j) / pivot;
					X.addVertices(i, rank, -quotient, j);
				}

				rank++;
			}
		}

		X.subMatrixRows(0, rank);
	}

	template <class R> CMatrix<R>* CHermiteNormalFormReducer<R>::reduceAndProduceUnitaryMatrix() {
		CMatrix<R>* U = new CMatrix<R>(CMatrix<R>::getUnitMatrix(X.getRowCount()));

		rank = 0;

		for (int j = 0; j < X.getColumnCount(); j++) {
			bool pivotFound = false;
			R pivot;
			for (int i = rank; i < X.getRowCount(); i++) {
				if (X(i, j) != 0) {
					if (!pivotFound) { // Found it now! :)
						X.swapVertices(rank, i);
						U->swapVertices(rank, i);
						pivotFound = true;
						pivot = X(rank, j);
						pivotColumns.push_back(j);

						// We must ensure it's positive
						if (pivot < 0) {
							X[rank] = -X[rank];
							(*U)[rank] = -(*U)[rank];
							pivot = -pivot;
						}

					} else {
						CExtendedExtendedGCDResults<R> xxGCD = extendedExtendedGCD(pivot, X(i, j));

						X.combineVertices(rank, i, xxGCD.leftBezoutFactor, xxGCD.rightBezoutFactor, xxGCD.leftLCMCofactor, -xxGCD.rightLCMCofactor);
						U->combineVertices(rank, i, xxGCD.leftBezoutFactor, xxGCD.rightBezoutFactor, xxGCD.leftLCMCofactor, -xxGCD.rightLCMCofactor);

						pivot = xxGCD.gcd;
					}
				}
			}

			// Now reduce the rows above the rank-row
			if (pivotFound) {
				for (int i = rank - 1; i >= 0; i--) {
					R quotient = X(i, j) / pivot;
					X.addVertices(i, rank, -quotient, j);
					U->addVertices(i, rank, -quotient, j);
				}

				rank++;
			}
		}

		X.subMatrixRows(0, rank);

		return U;
	}
}

#endif
