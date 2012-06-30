

#ifndef CVECTOR_H_
#define CVECTOR_H_

#include <vector>
using std::vector;

#include "../Structures.h"
#include <string>
using std::string;

#ifdef PPL_LINK
#include "ppl.h"
#endif

#include "../core/Vektor.h"

namespace AlgoTrans {
	template <class R> class CMatrix;

	template <class R> class CVector {
	protected:
		vector<R> elements;

	public:
		CVector() {};
		CVector(vector<R> iElements) : elements(iElements) {  }
		CVector(Vektor<R> v);

		void divideByGCDOfElements();
		CVector<R> getVectorDividedByGCDOfElements() const;

		void appendElement(R element) { elements.push_back(element); }

		CVector(int count, R iElements[]) {
			elements = vector<R>();
			for (int i = 0; i < count; i++) elements.push_back(iElements[i]);
		};

		CVector(int count) { elements = vector<R>(count); };
		CVector(const CVector<R>& iOriginal) { if (this != &iOriginal) *this = iOriginal; };

		R innerProduct(const CVector<R>& other) const;
		R normSquare() const;

		CVector<R>& operator = (const CVector<R>& other);
		CVector<R>& operator += (const CVector<R>& other);

		CVector<R>& operator << (const R& element);
		CVector<R> operator - () const;

		int getLength() const { return elements.size(); };

		// Operators
		R& operator [] (int element) { return elements[element]; }
		const R& operator [] (int element) const { return elements[element]; }

		operator CMatrix<R> () const; // Converts the vector to a row matrix
		operator Vektor<R> () const;

		CVector<R> getSubVector(int startElement, int length);

		static CVector<R> getUnitVector(int dim, int unitDim);
		static CVector<R> getZeroVector(int dim);
		static CVector<R> getConstantVector(int dim, R constant);

		string toString() const;
		string toStringHtml(bool horizontal = true) const;
		string toStringLatex(bool horizontal = false) const;
		void print() const;

		bool isZero() const;

#ifdef PPL_LINK
		operator Parma_Polyhedra_Library::Linear_Expression() const;
		static CVector<R> fromParmaConstraint(const Parma_Polyhedra_Library::Constraint& constraint);
		static CVector<R> fromParmaGenerator(const Parma_Polyhedra_Library::Generator& generator);
		static CVector<R> fromParmaRow(const Parma_Polyhedra_Library::Linear_Row& row);
#endif
	};

	template <class R>
	CVector<R>::CVector(Vektor<R> v) {
		for (int q = 0; q < v.getLength(); q++) elements.push_back(v[q]);
	}
}

#include "Vektor.cpp"

#endif

