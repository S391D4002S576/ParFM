#include "EchelonReduction.h"

namespace AlgoTrans {
	template <class R> void CEchelonReducer_VersionA<R>::reduce() {
		rank = 0;
		int m = X.generatingVertices.size();
		
		for (int j = 0; j < X.getSuperSpaceDimension(); j++) {
			// Search non-zero element in this column at row i (rank < i < m)
			bool found = false;
			for (int i0 = rank; (!found && (i0 < m)); i0++) {
				if (X(i0, j) != 0) found = true;
			}
			
			if (found) {
				pivotColumns.push_back(j);
				rank++;
				for (int i = m - 1; i >= rank; i--) {
					while (X(i, j) != 0) {
						R z = X(i - 1, j).getAbs() / X(i, j).getAbs();
						R s = sameSign(X(i, j), X(i - 1, j)) ? -z : z;
						X.addVertices(i - 1, i, s, j);
						X.swapVertices(i - 1, i);
					}
				}
			} 
		}
	}
	
	template <class R> void CEchelonReducer_VersionC<R>::reduce() {
		rank = 0;
		
		for (int j = 0; j < X.getSuperSpaceDimension(); j++) {
			for (int i = X.generatingVertices.size() - 1; i > rank; i--) {
				while (X(i, j) != 0) {
					R z = X(i - 1, j).getAbs() / X(i, j).getAbs();
					R s = sameSign(X(i, j), X(i - 1, j)) ? -z : z;
					X.addVertices(i - 1, i, s, j);
					X.swapVertices(i - 1, i);
				}
			}
			
			if (X(rank, j) != 0) {
				pivotColumns.push_back(j);
				rank++;
			}
		}
	}
	
	template <class R> void CEchelonReducer_VersionD<R>::reduce() {
		rank = 0;
		
		for (int j = 0; j < X.getSuperSpaceDimension(); j++) {
			bool pivotFound = false;
			R pivot;
			for (int i = rank; i < X.generatingVertices.size(); i++) {
				if (X(i, j) != 0) {
					if (!pivotFound) { // Found it now! :)
						X.swapVertices(rank, i);
						pivotFound = true;
						pivot = X(rank, j);
						pivotColumns.push_back(j);
					} else {
						CExtendedExtendedGCDResults<R> xxGCD = extendedExtendedGCD(pivot, X(i, j));

						printf("\n%s\n\n", X.toString().c_str());
						printf("[%s : %s : %s : %s : %s : %s]\n", 
								pivot.toString().c_str(), X(i, j).toString().c_str(), 
								xxGCD.leftBezoutFactor.toString().c_str(), xxGCD.rightBezoutFactor.toString().c_str(),
								xxGCD.leftLCMCofactor.toString().c_str(), xxGCD.rightLCMCofactor.toString().c_str());
						printf("[%i : %i]\n", 
								rank, i);
						X.combineVertices(rank, i, xxGCD.leftBezoutFactor, xxGCD.rightBezoutFactor,
								          xxGCD.leftLCMCofactor, -xxGCD.rightLCMCofactor);
						printf("\n%s\n\n", X.toString().c_str());
						
						pivot = xxGCD.gcd;
					}
				}
			}
			if (pivotFound) rank++;
		}
	}
	
	template <class R> void CScaledEchelonReducer<R>::reduce() {
		rank = 0;
		
		for (int j = 0; j < X.getSuperSpaceDimension(); j++) {
			bool pivotFound = false;
			R pivot;
			for (int i = rank; i < X.generatingVertices.size(); i++) {
				if (X(i, j) != 0) {
					if (!pivotFound) { // Found it now! :)
						X.swapVertices(rank, i);
						pivotFound = true;
						pivot = X(rank, j);
						pivotColumns.push_back(j);
					} else {
						CExtendedExtendedGCDResults<R> xxGCD = extendedExtendedGCD(pivot, X(i, j));

						printf("\n%s\n\n", X.toString().c_str());
						printf("[%s : %s : %s : %s : %s : %s]\n", 
								pivot.toString().c_str(), X(i, j).toString().c_str(), 
								xxGCD.leftBezoutFactor.toString().c_str(), xxGCD.rightBezoutFactor.toString().c_str(),
								xxGCD.leftLCMCofactor.toString().c_str(), xxGCD.rightLCMCofactor.toString().c_str());
						printf("[%i : %i]\n", 
								rank, i);
						X.combineVertices(rank, i, xxGCD.leftBezoutFactor, xxGCD.rightBezoutFactor,
								          xxGCD.leftLCMCofactor, -xxGCD.rightLCMCofactor);
						printf("\n%s\n\n", X.toString().c_str());
						
						pivot = xxGCD.gcd;
					}
				}
			}
			
			// Now reduce the rows above the rank-row
			if (pivotFound) {
				for (int i = rank - 1; i >= 0; i--) {
					CExtendedExtendedGCDResults<R> xxGCD = extendedExtendedGCD(pivot, X(i, j));

					printf("\n%s\n\n", X.toString().c_str());
					printf("[%s : %s : %s : %s : %s : %s]\n", 
							pivot.toString().c_str(), X(i, j).toString().c_str(), 
							xxGCD.leftBezoutFactor.toString().c_str(), xxGCD.rightBezoutFactor.toString().c_str(),
							xxGCD.leftLCMCofactor.toString().c_str(), xxGCD.rightLCMCofactor.toString().c_str());
					printf("[%i : %i]\n", 
							rank, i);
					X.combineVertex(i, rank, -xxGCD.rightLCMCofactor, xxGCD.leftLCMCofactor);
					printf("\n%s\n\n", X.toString().c_str());
				}
				
				rank++;
			}
		}
	}
}
