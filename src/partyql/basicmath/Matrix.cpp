#ifdef MATRIX_H_

#include "scalar/Integer.h"

namespace AlgoTrans {
#ifdef PPL_LINK
	typedef class Parma_Polyhedra_Library::Constraint_System ParmaConstraintSystem;
	typedef class Parma_Polyhedra_Library::Constraint ParmaConstraint;
	typedef class Parma_Polyhedra_Library::Linear_Expression ParmaLinearExpression;
#endif

	template <class R> CMatrix<R>::CMatrix(int iColumnCount, vector<CVector<R>*> iRows)
	: columnCount(iColumnCount),
	  rows(iRows)
	{
	}

	template <class R> void CMatrix<R>::swapVertices(int a, int b) {
		CVector<R>* v = this->rows[a];
		rows[a] = rows[b];
		rows[b] = v;
	}

	template <class R> void CMatrix<R>::addVertices(int a, int b, R factor) {
		addVertices(a, b, factor, 0);
	}

	template <class R> void CMatrix<R>::addVertices(int a, int b, R factor, int startColumn) {
		CMatrix& X = *this;
		for (int j = startColumn; j < getColumnCount(); j++) {
			X(a, j) += factor * X(b, j);
		}
	}

	template <class R> void CMatrix<R>::combineVertices(int a, int b, const R& aa, const R& ab, const R& ba, const R& bb, int startColumn) {
		CMatrix& X = *this;
		for (int j = startColumn; j < getColumnCount(); j++) {
			R rA = aa * X(a, j) + ab * X(b, j);
			R rB = ba * X(a, j) + bb * X(b, j);

			X(a, j) = rA;
			X(b, j) = rB;
		}
	}

	template <class R> void CMatrix<R>::combineVertices(int a, int b, const R& aa, const R& ab, const R& ba, const R& bb) {
		combineVertices(a, b, aa, ab, ba, bb, 0);
	}

	template <class R> void CMatrix<R>::combineVertex(int a, int b, const R& ca, const R& cb, int startColumn) {
		CMatrix& X = *this;
		for (int j = startColumn; j < getColumnCount(); j++) {
			X(a, j) = ca * X(a, j) + cb * X(b, j);
		}
	}

	template <class R> void CMatrix<R>::combineVertex(int a, int b, const R& ca, const R& cb) {
		combineVertex(a, b, ca, cb, 0);
	}

	template <class R> void CMatrix<R>::reduceToHermiteNormalForm() {
		CHermiteNormalFormReducer<R> reducer = CHermiteNormalFormReducer<R>(*this);

		reducer.reduce();
	}

	template <class R> CMatrix<R> CMatrix<R>::getHermiteNormalForm() const {
		CMatrix<R> result = CMatrix<R>(*this);

		result.reduceToHermiteNormalForm();

		return result;
	}

	template <class R> CMatrix<R> CMatrix<R>::getOrthogonalMatrix() const {
		CMatrix<R> X = CMatrix<R>(*this);

		X.transpose();

		CHermiteNormalFormReducer<R> hnfReducer = CHermiteNormalFormReducer<R>(X);
		CMatrix<R>* U = hnfReducer.reduceAndProduceUnitaryMatrix();

		U->subMatrixRows(hnfReducer.rank, U->getRowCount() - hnfReducer.rank);

		X = *U;

		delete U;

		return X;
	}

	template <class R> CMatrix<R> CMatrix<R>::getTransposedMatrix() const {
		CMatrix<R> result = *this;

		result.transpose();

		return result;
	}

	template <class R> int CMatrix<R>::getColumnCount() const {
		return columnCount;
	}

	template <class R> string CMatrix<R>::toString() const {
		int colCount = getColumnCount();
		if ((rows.size() == 0) || (colCount == 0)) {
			string result;

			result += "(Matrix has zero rows and/or columns)";

			return result;
		} else {
			int maxL = 0;
			vector<vector<string> > strings = vector<vector<string> >();
			for (unsigned int r = 0; r < rows.size(); r++) {
				strings.push_back(vector<string>(colCount));
				for (int c = 0; c < colCount; c++) {
					strings[r][c] = (*rows[r])[c].toString();
					int strL = strings[r][c].length();
					maxL = (strL > maxL) ? strL : maxL;
				}
			}

			string result;
			for (unsigned int r = 0; r < rows.size(); r++) {
				result += "[ ";
				for (int c = 0; c < colCount; c++) {
					int spaceCount = (maxL - strings[r][c].length()) - (c == 0 ? 1 : 0);
					for (int s = spaceCount; s >= 0; s--) result += " ";
					result += strings[r][c];
				}
				result += " ]\n";
			}

			return result;
		}
	}

	template <class R> string CMatrix<R>::toStringHtml() const {
		int colCount = getColumnCount();
		if ((rows.size() == 0) || (colCount == 0)) {
			string result;

			result += "(Matrix has zero rows and/or columns)";

			return result;
		} else {
			string result;
			result += "<table cellspacing=\"0\" border=\"1\">";
			for (unsigned int r = 0; r < rows.size(); r++) {
				result += "<tr>";
				for (int c = 0; c < colCount; c++) {
					result += "<td><center>";
					result += (*rows[r])[c].toString();
					result += "</center></td>";
				}
				result += "</tr>";
			}
			result += "</table>";

			return result;
		}
	}

	template <class R> string CMatrix<R>::toStringLatex() const {
		int colCount = getColumnCount();
		if ((rows.size() == 0) || (colCount == 0)) {
			string result;

			result += "(Matrix has zero rows and/or columns)";

			return result;
		} else {
			string result;
			result += "\\left[ \\begin{array}{";
			for (int c = 0; c < colCount; c++) result += "c";
			result += "}\n";
			for (unsigned int r = 0; r < rows.size(); r++) {
				if (r != 0) result += "\\\\\n";
				for (int c = 0; c < colCount; c++) {
					if (c != 0) result += "&";
					result += (*rows[r])[c].toString();
				}
			}
			result += "\\end{array} \\right]";

			return result;
		}
	}

	template <class R> void CMatrix<R>::transpose() {
		vector<CVector<R>*> newVertices = vector<CVector<R>*>();

		for (int c = 0; c < getColumnCount(); c++) {
			newVertices.push_back(new CVector<R>(rows.size()));
		}

		for (unsigned int r = 0; r < rows.size(); r++) {
			for (int c = 0; c < getColumnCount(); c++) {
				(*newVertices[c])[r] = (*rows[r])[c];
			}
		}

		for (unsigned int r = 0; r < rows.size(); r++) {
			delete rows[r];
		}

		columnCount = rows.size();
		rows = newVertices;
	}

	template <class R> void CMatrix<R>::subMatrixRows(int startRow, int count) {
		vector<CVector<R>*> newVertices = vector<CVector<R>*>(count);

		for (int q = startRow; q < startRow + count; q++) {
			newVertices[q - startRow] = rows[q];
		}

		// Delete discarded vertices
		for (int q = 0; q < startRow; q++) delete rows[q];
		for (unsigned int q = startRow + count; q < rows.size(); q++) delete rows[q];

		rows = newVertices;
	}

	template <class R> CMatrix<R> CMatrix<R>::getSubMatrix(int startRow, int rowCount, int startCol, int colCount) const {
		ASSERT((0 <= startRow) && (startRow + rowCount <= getRowCount()));
		ASSERT((0 <= startCol) && (startCol + colCount <= getColumnCount()));

		CMatrix<R> result = CMatrix<R>(colCount);

		for (int q = startRow; q < startRow + rowCount; q++) {
			result.addRow(rows[q]->getSubVector(startCol, colCount));
		}

		return result;
	}

	template <class R> CMatrix<R> CMatrix<R>::getSubMatrixRows(int startRow, int count) const {
		return getSubMatrix(startRow, count, 0, getColumnCount());
	}

	template <class R> CMatrix<R> CMatrix<R>::getSubMatrixColumns(int startColumn, int count) const {
		return getSubMatrix(0, getRowCount(), startColumn, count);
	}

	template <class R> bool operator == (const CMatrix<R>& a, const CMatrix<R>& b) {
		bool equal = (a.getColumnCount() == b.getColumnCount());

		for (unsigned int q = 0; (q < a.rows.size()) && equal; q++) {
			equal = equal && (a[q] == b[q]);
		}

		return equal;
	}

	template <class R> CMatrix<R> operator - (const CMatrix<R>& a, const CMatrix<R>& b) {
		ASSERT((a.getColumnCount() == b.getColumnCount()) && (a.getRowCount() == b.getRowCount()));

		CMatrix<R> result = CMatrix<R>(a.getColumnCount());

		for (int q = 0; q < a.getRowCount(); q++) {
			result.addRow(a[q] - b[q]);
		}

		return result;
	}

	template <class R> CMatrix<R> operator + (const CMatrix<R>& a, const CMatrix<R>& b) {
		ASSERT((a.getColumnCount() == b.getColumnCount()) && (a.getRowCount() == b.getRowCount()));

		CMatrix<R> result = CMatrix<R>(a.getColumnCount());

		for (int q = 0; q < a.getRowCount(); q++) {
			result.addRow(a[q] + b[q]);
		}

		return result;
	}

	template <class R> CMatrix<R> operator * (const CMatrix<R>& a, const CMatrix<R>& b) {
		ASSERT(a.getColumnCount() == b.getRowCount());

		CMatrix<R> result = CMatrix<R>(b.getColumnCount());
		for (int q = 0; q < a.getRowCount(); q++) {
			CVector<R> row = CVector<R>();
			for (int r = 0; r < b.getColumnCount(); r++) {
				R sum = 0;
				for (int k = 0; k < a.getColumnCount(); k++) sum += a[q][k] * b[k][r];
				row.appendElement(sum);
			}
			result.addRow(row);
		}

		return result;
	}

	// Null matrix
	template <class R> CMatrix<R> operator ! (const CMatrix<R>& a) {
		return a.getOrthogonalMatrix();
	}

	template <class R> CMatrix<R>& CMatrix<R>::operator = (const CMatrix<R>& other) {
		if (this != &other) {
			for (unsigned int q = 0; q < rows.size(); q++) delete rows[q];

			columnCount = other.columnCount;
			rows = vector<CVector<R>*>();
			for (int i = 0; i < other.getRowCount(); i++) {
				rows.push_back(new CVector<R>(other[i]));
			}
		}

		return *this;
	}

	template <class R> void CMatrix<R>::print() const {
		printf("%s\n", toString().c_str());
	}

	/*template <class R> void CMatrix<R>::addRow(CVector<R>* v) {
		ASSERT(v->getLength() == getColumnCount());

		rows.push_back(v);
	}*/

	template <class R> void CMatrix<R>::addRow(const CVector<R>& v) {
		ASSERT_EQ(CInteger(v.getLength()), CInteger(getColumnCount()));

		rows.push_back(new CVector<R>(v));
	}

	template <class R> void CMatrix<R>::addRow(CVector<R>* v) {
		ASSERT(v->getLength() == getColumnCount());

		rows.push_back(v);
	}

	/*template <class R> CMatrix<R>* CMatrix<R>::getUnitMatrix(int dim) {
		CMatrix<R>* result = new CMatrix<R>(dim);

		for (int q = 0; q < dim; q++) {
			result->addRow(new CVector<R>(CVector<R>::getUnitVector(dim, q)));
		}

		return result;
	}*/

	template <class R> CMatrix<R> operator << (const CMatrix<R>& a, const CMatrix<R>& b) {
		if ((a.getColumnCount() == 0) || (a.getRowCount() == 0)) return b;
		if ((b.getColumnCount() == 0) || (b.getRowCount() == 0)) return a;

		ASSERT(a.getRowCount() == b.getRowCount());
		CMatrix<R> result = CMatrix<R>(a.getColumnCount() + b.getColumnCount());

		for (int r = 0; r < a.getRowCount(); r++) {
			result.addRow(a[r] << b[r]);
		}

		return result;
	}

	/*template <class R> CVector<R> operator << (const CMatrix<R>& a, const CVector<R>& b) {
		ASSERT(a.getRowCount() == 1);

		return a[0] << b;
	}

	template <class R> CVector<R> operator << (const CVector<R>& a, const CMatrix<R>& b) {
		ASSERT(a.getRowCount() == 1);

		return a << b[0];
	}*/

	template <class R> CMatrix<R>& CMatrix<R>::operator >>= (const CMatrix<R>& m) {
		for (int r = 0; r < m.getRowCount(); r++) {
			ASSERT(m[r].getLength() == getColumnCount());

			addRow(CVector<R>(m[r]));
		}

		return *this;
	}

	template <class R> CMatrix<R>& CMatrix<R>::operator >>= (const CVector<R>& v) {
		ASSERT(v.getLength() == getColumnCount());

		addRow(v);

		return *this;
	}

	template <class R> CMatrix<R> operator >> (const CMatrix<R>& a, const CMatrix<R>& b) {
		if (a.getRowCount() == 0) return b;
		if (b.getRowCount() == 0) return a;
		ASSERT_EQ(CInteger(a.getColumnCount()), CInteger(b.getColumnCount()));

		CMatrix<R> result = CMatrix<R>(a.getColumnCount());

		for (int r = 0; r < a.getRowCount(); r++) result.addRow(a[r]);
		for (int r = 0; r < b.getRowCount(); r++) result.addRow(b[r]);

		return result;
	}

	template <class R> CMatrix<R> operator >> (const CMatrix<R>& a, const CVector<R>& b) {
		ASSERT_EQ(CInteger(a.getColumnCount()), CInteger(b.getLength()));
		CMatrix<R> result = CMatrix<R>(a.getColumnCount());

		for (int r = 0; r < a.getRowCount(); r++) result.addRow(a[r]);
		result.addRow(b);

		return result;
	}

	template <class R> CMatrix<R> operator >> (const CVector<R>& a, const CMatrix<R>& b) {
		ASSERT_EQ(CInteger(a.getLength()), CInteger(b.getColumnCount()));
		CMatrix<R> result = CMatrix<R>(a.getLength());

		result.addRow(new CVector<R>(a));
		for (int r = 0; r < b.getRowCount(); r++) result.addRow(b[r]);

		return result;
	}

	template <class R> CMatrix<R> CMatrix<R>::getZeroMatrix(int rows, int cols) {
		CMatrix<R> result = CMatrix<R>(cols);

		for (int q = 0; q < rows; q++) result.addRow(CVector<R>::getZeroVector(cols));

		return result;
	}

	template <class R> CMatrix<R> CMatrix<R>::getUnitMatrix(int dim) {
		CMatrix<R> result = CMatrix<R>(dim);

		for (int q = 0; q < dim; q++) result.addRow(CVector<R>::getUnitVector(dim, q));

		return result;
	}

	template <class R> CMatrix<R> operator - (const CMatrix<R> a) {
		CMatrix<R> result = CMatrix<R>(a.getColumnCount());

		for (int q = 0; q < a.getRowCount(); q++) result.addRow(-a[q]);

		return result;
	}

	template <class R>
	CMatrix<R>::~CMatrix() {
		for (unsigned int q = 0; q < rows.size(); q++) delete rows[q];
	}

	template <class R>
	CMatrix<R> CMatrix<R>::polyhedrallyDecomposeConstraintMatrix(CLattice<R> core) const {
		const CMatrix<R>& X = *this;

		CMatrix<R> coreMatrix = CMatrix<R>(core.getMatrix());
		CHermiteNormalFormReducer<R> reducer = CHermiteNormalFormReducer<R>(coreMatrix);
		reducer.reduce();

		vector<int> usefulPivots;
		vector<R> multipliers;
		for (unsigned int p = 0; p < reducer.pivotColumns.size(); p++) {
			if (!coreMatrix[p][reducer.pivotColumns[p]].isOne()) {
				usefulPivots.push_back(reducer.pivotColumns[p]);
				multipliers.push_back(coreMatrix[p][reducer.pivotColumns[p]]);
			}
		}

		CMatrix<R> newConstraintMatrix = CMatrix<R>(getColumnCount() + multipliers.size());
		for (int r = 0; r < getRowCount(); r++) {
			unsigned int nextPivot = 0;
			CVector<R> row;
			for (int c = 0; c < getColumnCount(); c++) {
				int nextPivIx = (nextPivot < usefulPivots.size()) ? usefulPivots[nextPivot] : -1;
				if (nextPivIx == c) {
					row.appendElement(multipliers[nextPivIx] * X(r, c));
					row.appendElement(X(r, c));
					nextPivot++;
				} else {
					row.appendElement(X(r, c));
				}
			}
			newConstraintMatrix.addRow(row);
		}

		return newConstraintMatrix;
	}

#ifdef PPL_LINK
	template <class R>
	CMatrix<R>::operator ParmaConstraintSystem() const {
		ParmaConstraintSystem result = ParmaConstraintSystem();

		for (int q = 0; q < getRowCount(); q++) {
			result.insert(((ParmaLinearExpression) (*this)[q]) >= 0);
		}

		return result;
	}

	typedef class Parma_Polyhedra_Library::Constraint_System ParmaConstraintSystem;
	typedef class Parma_Polyhedra_Library::Generator_System ParmaGeneratorSystem;
	typedef class Parma_Polyhedra_Library::Linear_System ParmaLinearSystem;
	typedef ParmaConstraintSystem::const_iterator ConstraintCIterator;
	typedef ParmaGeneratorSystem::const_iterator GeneratorCIterator;
	typedef ParmaLinearSystem::const_iterator LinearCIterator;

	template <class R>
	CMatrix<R> CMatrix<R>::fromParmaConstraintSystem(const ParmaConstraintSystem & constraintSystem) {
		CMatrix<R> result = CMatrix<R>(constraintSystem.space_dimension() + 1);

		for (ConstraintCIterator q = constraintSystem.begin(); q != constraintSystem.end(); ++q) {
			CVector<R> v = CVector<R>::fromParmaConstraint(*q);
			result.addRow(v);
			if ((*q).type() == Parma_Polyhedra_Library::Constraint::EQUALITY) result.addRow(-v);
		}

		return result;
	}

	template <class R>
	CMatrix<R> CMatrix<R>::fromParmaGeneratorSystem(const Parma_Polyhedra_Library::Generator_System & generatorSystem) {
		CMatrix<R> result = CMatrix<R>(generatorSystem.space_dimension() + 1);

		//generatorSystem.ascii_dump(std::cout);
		for (GeneratorCIterator q = generatorSystem.begin(); q != generatorSystem.end(); ++q) {
			CVector<R> v = CVector<R>::fromParmaGenerator(*q);
			result.addRow(v);
			if ((*q).type() == Parma_Polyhedra_Library::Generator::LINE) result.addRow(-v);
		}

		return result;
	}
#endif
}

#endif /* MODULE_CPP_ */
