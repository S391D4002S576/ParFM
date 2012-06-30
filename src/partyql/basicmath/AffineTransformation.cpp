#ifndef AFFINETRANSFORMATION_CPP_
#define AFFINETRANSFORMATION_CPP_

#include "AffineTransformation.h"

namespace AlgoTrans {
	template <class R> CAffineTransformation<R>::CAffineTransformation(int iSourceSpaceDimension, vector<CVector<R>*> iTransformationVectors)
	: CMatrix<R>(iSourceSpaceDimension + 1, iTransformationVectors)
	{
	}

	template <class R> string CAffineTransformation<R>::toStringLatex(vector<string>* dimensionNames) const {
		if (dimensionNames == NULL) {
			return CMatrix<R>::toStringLatex();
		} else {
			if (CMatrix<R>::getRowCount() == 0) {
				return "(Affine transformation to 0-dimensional space)";
			} else {
				const CMatrix<R>& X = *this;
				string result;
				if (CMatrix<R>::getRowCount() > 1) result += "\\left[ \\begin{array}{c}\n";
				int colCount = CMatrix<R>::getColumnCount();
				for ( int r = 0; r < CMatrix<R>::getRowCount(); r++) {
					if (r != 0) result += "\\\\\n";
					bool firstVal = true;
					for (int c = 0; c < colCount; c++) { R v = X(r, c);
						if (v != 0) {
							string dimName = (c < (int) dimensionNames->size()) ? (*dimensionNames)[c] : "";
							R a = v.getAbs();
							result += (v < 0) ? (firstVal ? "-" : " - ")
									          : (firstVal ? "" : " + ");
							if ((a != 1) || (dimName == "")) result += a.toString();
							result += dimName;
							firstVal = false;
						}
					}
					if (firstVal) result += "0";
				}
				if (CMatrix<R>::getRowCount() > 1) result += "\\end{array} \\right]";

				return result;
			}
		}
	}

	template <class R> string CAffineTransformation<R>::toStringHtml(vector<string>* dimensionNames) const {
		if (dimensionNames == NULL) {
			return CMatrix<R>::toStringHtml();
		} else {
			if (CMatrix<R>::getRowCount() == 0) {
				return "(Affine transformation to 0-dimensional space)";
			} else {
				const CMatrix<R>& X = *this;
				string result;
				if (CMatrix<R>::getRowCount() > 1) result += "<table><tr><td>";
				int colCount = CMatrix<R>::getColumnCount();
				for ( int r = 0; r < CMatrix<R>::getRowCount(); r++) {
					if (r != 0) result += "</td></tr><tr><td>";
					bool firstVal = true;
					for (int c = 0; c < colCount; c++) { R v = X(r, c);
						if (v != (R) 0) {
							string dimName = (c < (int) dimensionNames->size()) ? (*dimensionNames)[c] : "";
							R a = v.getAbs();
							result += (v < (R) 0) ? (firstVal ? "-" : " - ")
									          : (firstVal ? "" : " + ");
							if ((a != (R) 1) || (dimName == "")) result += a.toString();
							result += dimName;
							firstVal = false;
						}
					}
					if (firstVal) result += "0";
				}
				if (CMatrix<R>::getRowCount() > 1) result += "</td></tr></table>";

				return result;
			}
		}
	}

	template <class R> bool operator == (const CAffineTransformation<R> a, const CAffineTransformation<R> b) {
		return a.matrix == b.matrix;
	}

	template <class R> CAffineTransformation<R>& CAffineTransformation<R>::operator = (const CAffineTransformation<R>& other) {
		if (this != &other) {
			CMatrix<R>::operator=(other);
		}

		return *this;
	}

	template <class R> void CAffineTransformation<R>::addTransformationRow(CVector<R>* v) {
		CMatrix<R>::addRow(*v); delete v;
	}

	template <class R>
	CAffineTransformation<R> CAffineTransformation<R>::polyhedralLatticeDecomposition(const CLattice<R>& core) {
		return CAffineTransformation<R>(CMatrix<R>::polyhedrallyDecomposeConstraintMatrix(core));
	}
}

#endif /* AFFINETRANSFORMATION_CPP_ */
