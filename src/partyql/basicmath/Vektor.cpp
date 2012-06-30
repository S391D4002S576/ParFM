

#ifndef VEKTOR_CPP_
#define VEKTOR_CPP_

#include "Vektor.h"
#include "../cute/cute.h"

#ifdef PPL_LINK
#include "ppl.h"
#endif

namespace AlgoTrans {
#ifdef PPL_LINK
	typedef class Parma_Polyhedra_Library::Linear_Expression ParmaLinearExpression;
	typedef class Parma_Polyhedra_Library::Variable ParmaVariable;
#endif

	template <class R> bool CVector<R>::isZero() const {
		bool result = true;

		for (int q = 0; q < getLength(); q++) result = result && elements[q].isZero();

		return result;
	}

	template <class R> string CVector<R>::toString() const {
		int maxL = 0;
		string strings[getLength()];
		for (int q = 0; q < getLength(); q++) {
			strings[q] = elements[q].toString();
			int strL = strings[q].length();
			maxL = (strL > maxL) ? strL : maxL;
		}

		string result;
		for (int q = 0; q < getLength(); q++) {
			for (int s = (maxL - strings[q].length()); s >= 0; s--) result += " ";
			result += strings[q];
		}

		return result;
	}

	template <class R> string CVector<R>::toStringLatex(bool horizontal) const {
		string result;
		result += "\\left[ \\begin{array}{";
		if (horizontal) {
			for (int c = 0; c < getLength(); c++) result += "c";
			result += "}\n";
			for (int q = 0; q < getLength(); q++) {
				if (q != 0) result += "&";
				result += elements[q].toString();
			}
		} else {
			result += "c}\n";
			for (int q = 0; q < getLength(); q++) {
				if (q != 0) result += "\\\\";
				result += elements[q].toString();
			}
		}
		result += "\\end{array} \\right] ";

		return result;
	}

	template <class R> CVector<R>& CVector<R>::operator = (const CVector<R>& other) {
		if (this != &other) {
			elements = vector<R>();
			for (unsigned int i = 0; i < other.elements.size(); i++) {
				elements.push_back(other.elements[i]);
			}
		}

		return *this;
	}

	template <class R> CVector<R>::operator CMatrix<R> () const { // Converts the vector to a row matrix
		CMatrix<R> result = CMatrix<R>(getLength());

		result.addRow(*this);

		return result;
	}

	template <class R> CVector<R> operator + (const CVector<R>& a, const CVector<R>& b) {
		ASSERT(a.getLength() == b.getLength());

		vector<R> result = vector<R>();

		for (int q = 0; q < a.getLength(); q++) result.push_back(a[q] + b[q]);

		return result;
	};

	template <class R> CVector<R> operator * (const R& f, const CVector<R>& v) {
		vector<R> result = vector<R>();

		for (int q = 0; q < v.getLength(); q++) result.push_back(f * v[q]);

		return result;
	};

	template <class R> CVector<R> operator - (const CVector<R>& a, const CVector<R>& b) {
		ASSERT(a.getLength() == b.getLength());

		CVector<R> result = CVector<R>(a.getLength());

		for (int q = 0; q < a.getLength(); q++) result[q] = (a[q] - b[q]);

		return result;
	};

	template <class R> CVector<R> CVector<R>::operator - () const {
		CVector<R> result = CVector<R>(getLength());

		for (int q = 0; q < getLength(); q++) result[q] = -elements[q];

		return result;
	};

	template <class R> CVector<R> operator << (const CVector<R>& a, const CVector<R>& b) {
		CVector<R> result = CVector<R>();

		for (int q = 0; q < a.getLength(); q++) result.appendElement(a[q]);
		for (int q = 0; q < b.getLength(); q++) result.appendElement(b[q]);

		return result;
	};

	template <class R> bool operator == (const CVector<R>& a, const CVector<R>& b) {
		bool equal = (a.getLength() == b.getLength());

		for (int e = 0; equal && (e < a.getLength()); e++) {
			equal = equal && (a[e] == b[e]);
		}

		return equal;
	};

	template <class R> bool operator != (const CVector<R>& a, const CVector<R>& b) {
		return !(a == b);
	}

	template <class R> CVector<R> operator - (const CVector<R>& a) {
		vector<R> result = vector<R>();

		for (int q = 0; q < a.getLength(); q++) result.push_back(-a[q]);

		return result;
	}

	template <class R> CVector<R> CVector<R>::getUnitVector(int dim, int unitDim) {
		CVector<R> result = CVector<R>();

		for (int q = 0; q < dim; q++) {
			result.elements.push_back((q == unitDim) ? 1 : 0);
		}

		return result;
	}

	template <class R> CVector<R> CVector<R>::getZeroVector(int dim) {
		return getConstantVector(dim, 0);
	}

	template <class R> CVector<R> CVector<R>::getConstantVector(int dim, R constant) {
		CVector<R> result = CVector<R>(dim);

		for (int q = 0; q < dim; q++) result[q] = constant;

		return result;
	}

	template <class R> void CVector<R>::print() const {
		printf("%s\n", toString().c_str());
	}

	template <class R> CVector<R>& CVector<R>::operator << (const R& element) {
		appendElement(element);

		return *this;
	}

	template <class G> CVector<G> operator << (const G& a, const CVector<G>& b) { // hor concatenation
		CVector<G> result;

		result.appendElement(a);
		for (int q = 0; q < b.getLength(); q++) result.appendElement(b[q]);

		return result;
	}

	template <class G> bool operator < (const CVector<G>& a, const CVector<G>& b) { // lexicographical ordering
		for (int q = 0; q < b.getLength(); q++) {
			if (a[q] < b[q]) {
				return true;
			} else if (a[q] > b[q]) {
				return false;
			}
		}

		return false;
	}

	template <class R> CVector<R> CVector<R>::getSubVector(int startElement, int length) {
		ASSERT((0 <= startElement) && (startElement + length <= getLength()));

		CVector<R> result = CVector<R>();

		for (int q = startElement; q < startElement + length; q++) {
			result.elements.push_back((*this)[q]);
		}

		return result;
	}

	template <class R> void CVector<R>::divideByGCDOfElements() {
		R g = elements[0];
		for (int c = 1; c < getLength(); c++) g = gcd(g, elements[c]);

		if (g != 0) {
			g = (g > 0) ? g : -g;
			for (int c = 0; c < getLength(); c++) elements[c] = elements[c] / g;
		}
	}

	template <class R> CVector<R> CVector<R>::getVectorDividedByGCDOfElements() const {
		CVector<R> result = *this;

		result.divideByGCDOfElements();

		return result;
	}

	/*template <class R>
	CVector<R>::operator ParmaLinearExpression() const {
		ParmaLinearExpression result = ParmaLinearExpression();

		for (unsigned int q = 0; q < elements.size(); q++) {
			if (q > 0) {
				result += elements[q].getMPZClass() * ParmaVariable(q - 1);
			} else {
				result += elements[q].getMPZClass();
			}
		}

		return result;
	}

	template <class R>
	CVector<R> CVector<R>::fromParmaConstraint(const Parma_Polyhedra_Library::Constraint& constraint) {
		CVector<R> result = CVector<R>(constraint.space_dimension() + 1);

		for (int q = 1; q < result.getLength(); q++) {
			result[q] = R(constraint.coefficient(Parma_Polyhedra_Library::Variable(q - 1)));
		}
		result[0] = R(constraint.inhomogeneous_term());

		return result;
	}

	template <class R>
	CVector<R> CVector<R>::fromParmaRow(const Parma_Polyhedra_Library::Linear_Row& row) {
		CVector<R> result = CVector<R>(row.space_dimension() + 1);

		for (int q = 1; q < result.getLength(); q++) {
			result[q] = R(row.coefficient(q - 1));
		}
		result[0] = R(row.inhomogeneous_term());

		return result;
	}

	template <class R>
	CVector<R> CVector<R>::fromParmaGenerator(const Parma_Polyhedra_Library::Generator& generator) {
		CVector<R> result = CVector<R>(generator.space_dimension() + 1);

		for (unsigned int q = 0; q < generator.space_dimension(); q++) {
			result[1 + q] = generator.coefficient(Parma_Polyhedra_Library::Variable(q));
		}
		if (generator.is_point() || generator.is_closure_point()) {
			result[0] = 1;
		} else if (generator.is_ray()) {
			result[0]  = 0;
		} else {
			result[0] = generator.divisor();
		}

		return result;
	}*/

	template <class R>
	R CVector<R>::innerProduct(const CVector<R>& other) const {
		ASSERT(getLength() == other.getLength());

		R result = 0;

		for (int q = 0; q < other.getLength(); q++) result += (*this)[q] * other[q];

		return result;
	}

	template <class R>
	CVector<R>& CVector<R>::operator += (const CVector<R>& other) {
		for (int q = 0; q < getLength(); q++) elements[q] += other[q];

		return *this;
	}

	template <class R>
	R CVector<R>::normSquare() const { return innerProduct(*this); }

	template <class R>
	std::string CVector<R>::toStringHtml(bool horizontal) const {
		string result = "<table><tr><td>Vec</td><td>";
		if ((elements.size() != 0)) {
			result += "<table cellspacing=\"0\" border=\"1\">";
			if (horizontal) result += "<tr>";
			for (unsigned int r = 0; r < elements.size(); r++) {
				if (!horizontal) result += "<tr>";
				result += "<td>" + elements[r].toString() + "</td>";
				if (!horizontal) result += "</tr>";
			}
			if (horizontal) result += "</tr>";
			result += "</table>";
		}

		return result + "</td></tr></table>";
	}

	template <class R>
	CVector<R>::operator Vektor<R> () const {
		Vektor<R> result = Vektor<R>(getLength());

		for (int q = 0; q < getLength(); q++) result[q] = elements[q];

		return result;
	}
}

#endif /*VEKTOR_CPP_*/


