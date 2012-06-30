#ifndef AFFINETRANSFORMATION_H_
#define AFFINETRANSFORMATION_H_

#include <vector>
using std::vector;

#include <string>
using std::string;
#include <list>
using std::list;

#include <algorithm>
using std::min;

#include "Vektor.h"
#include "Matrix.h"

namespace AlgoTrans {
	template <class R> class CAffineTransformation: public CMatrix<R> {
		template <class G> friend bool operator == (const CAffineTransformation<G> a, const CAffineTransformation<G> b);
	private:
	protected:
	public:
		const CMatrix<R>& getMatrix() const { return *this; }
		
		CAffineTransformation() { }
		CAffineTransformation(int iSourceSpaceDimension) : CMatrix<R>(iSourceSpaceDimension + 1) { }
		CAffineTransformation(int iSourceSpaceDimension, vector<CVector<R>*> iTransformationVectors);
		CAffineTransformation(const CMatrix<R>& iMatrix) : CMatrix<R>(iMatrix) { }
		CAffineTransformation(const CAffineTransformation<R>& iOriginal) : CMatrix<R>(iOriginal.getSourceSpaceDimension()) { if (this != &iOriginal) *this = iOriginal; };
		
		void addTransformationRow(CVector<R>* v);
		
		int getSourceSpaceDimension() const { return CMatrix<R>::getColumnCount() - 1; };
		int getDestinationSpaceDimension() const { return CMatrix<R>::getRowCount(); };
		
		string toStringLatex(vector<string>* dimensionNames = NULL) const;
		string toStringHtml(vector<string>* dimensionNames = NULL) const;
		
		CAffineTransformation<R>& operator = (const CAffineTransformation<R>& other);
		
		CAffineTransformation<R> polyhedralLatticeDecomposition(const CLattice<R>& core);
	};
}

#include "AffineTransformation.cpp"

#endif /* AFFINETRANSFORMATION_H_ */
