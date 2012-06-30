#ifndef MATRIX_H_
#define MATRIX_H_

#include <vector>
using std::vector;

#include <list>
using std::list;

#include <algorithm>
using std::min;

#include "Vektor.h"
#include "factorization/HermiteNormalFormReduction.h"

#include "../cute/cute.h"

//#include "ppl.h"
#include "VectorDescriptor.h"

namespace AlgoTrans {
	template <class R> class CLattice;

	template <class R> class CMatrix {
		//template <class G> friend CMatrix<G> operator (const vector<CVector<G>* >);
		template <class G> friend class CHermiteNormalFormReducer;
		template <class G> friend bool operator == (const CMatrix<G>& a, const CMatrix<G>& b);
		template <class G> friend CMatrix<G> operator - (const CMatrix<G>& a, const CMatrix<G>& b);
		template <class G> friend CMatrix<G> operator + (const CMatrix<G>& a, const CMatrix<G>& b);
		template <class G> friend CMatrix<G> operator * (const CMatrix<G>& a, const CMatrix<G>& b);
		template <class G> friend CMatrix<G> operator ! (const CMatrix<G>& a); // Unitary nullspace
		template <class G> friend CMatrix<G> operator << (const CMatrix<G>& a, const CMatrix<G>& b); // hor concatenation
		//template <class G> friend CMatrix<G> operator << (const CMatrix<G>& a, const CVector<G>& b); // hor concatenation
		//template <class G> friend CMatrix<G> operator << (const CVector<G>& a, const CMatrix<G>& b); // hor concatenation

		template <class G> friend CMatrix<G> operator >> (const CMatrix<G>& a, const CMatrix<G>& b); // ver concatenation
		template <class G> friend CMatrix<G> operator >> (const CMatrix<G>& a, const CVector<G>& b); // ver concatenation
		template <class G> friend CMatrix<G> operator >> (const CVector<G>& a, const CMatrix<G>& b); // ver concatenation
	private:
		void swapVertices(int a, int b);

		// Replace vector <a> with <a> + factor * <b>
		void addVertices(int a, int b, R factor);
		void addVertices(int a, int b, R factor, int startColumn);

		// Replace vector <a> with aa * <a> + ab * <b>
		//         vector <b> with ba * <a> + bb * <b>
		void combineVertices(int a, int b, const R& aa, const R& ab, const R& ba, const R& bb);
		void combineVertices(int a, int b, const R& aa, const R& ab, const R& ba, const R& bb, int startColumn);

		// Replace vector <a> with ca * <a> + cb * <b>
		void combineVertex(int a, int b, const R& ca, const R& cb, int startColumn);
		void combineVertex(int a, int b, const R& ca, const R& cb);

	protected:

		int columnCount; // We need to store this value since otherwise we don't know the dimension of the orthogonal module
		vector<CVector<R>*> rows;
		R& operator () (int i, int j) { return (*rows[i])[j]; };

	public:
		~CMatrix<R>();

		const R& operator () (int i, int j) const { return (*rows[i])[j]; };

		CVector<R>& operator [] (int v) { return *rows[v]; };
		CVector<R>& operator [] (int v) const { return *rows[v]; };

		CMatrix<R>& operator >>= (const CMatrix<R>& m);
		CMatrix<R>& operator >>= (const CVector<R>& v);

		CMatrix() : columnCount(-1) { }
		CMatrix(int iColumnCount) : columnCount(iColumnCount) { }
		CMatrix(int iColumnCount, vector<CVector<R>*> iRows);
		CMatrix(const CMatrix<R>& iOriginal) { if (this != &iOriginal) *this = iOriginal; };

		//void addRow(CVector<R>* v) { addRow(*v); delete v; };
		void addRow(const CVector<R>& v);
		void addRow(CVector<R>* v);

		void subMatrixRows(int startRow, int count);
		CMatrix<R> getSubMatrix(int startRow, int rowCount, int startCol, int colCount) const;
		CMatrix<R> getSubMatrixRows(int startRow, int count) const;
		CMatrix<R> getSubMatrixColumns(int startColumn, int count) const;
		CMatrix<R> getColumn(int column) const { return getSubMatrixColumns(column, 1); }

		void reduceToHermiteNormalForm();
		CMatrix<R> getHermiteNormalForm() const;

		CMatrix<R> getOrthogonalMatrix() const;

		void transpose();
		CMatrix<R> getTransposedMatrix() const;

		int getColumnCount() const;
		int getRowCount() const { return rows.size(); };

		string toString() const;
		void print() const;
		string toStringHtml() const;
		string toStringLatex() const;

		CMatrix<R>& operator = (const CMatrix<R>& other);

		static CMatrix<R> getUnitMatrix(int dim);
		static CMatrix<R> getZeroMatrix(int rows, int cols);

		CMatrix polyhedrallyDecomposeConstraintMatrix(CLattice<R> core) const;

#ifdef PPL_LINK
		operator Parma_Polyhedra_Library::Constraint_System() const;
		static CMatrix<R> fromParmaConstraintSystem(const Parma_Polyhedra_Library::Constraint_System & constraintSystem);
		static CMatrix<R> fromParmaGeneratorSystem(const Parma_Polyhedra_Library::Generator_System & generatorSystem);
#endif
	};
}

#include "Matrix.cpp"

#endif /* MATRIX_H_ */
