#ifndef CVEKTOR_H_
#define CVEKTOR_H_

#include <vector>
using std::vector;

#include <string>
using std::string;

namespace AlgoTrans {
	template <class R> class Vektor {
	protected:
		vector<R> elements;

	public:
		Vektor() {};
		Vektor(vector<R> iElements) : elements(iElements) {  }

		void appendElement(R element) { elements.push_back(element); }

		Vektor(int count, R iElements[]) { elements = vector<R>(); for (int i = 0; i < count; i++) elements.push_back(iElements[i]); };

		Vektor(int count) { elements = vector<R>(count); };
		Vektor(int count, R defaultValue) { elements = vector<R>(count, defaultValue); };

		Vektor(const Vektor<R>& iOriginal) { if (this != &iOriginal) *this = iOriginal; };

		Vektor<R>& operator = (const Vektor<R>& other) {
			if (this != &other) {
				elements = other.elements;
			}

			return *this;
		}

		Vektor<R>& operator += (const Vektor<R>& other);

		Vektor<R>& operator << (const R& element);

		Vektor<R> operator -() const;

		int getLength() const { return elements.size(); };

		R& operator [] (int element) { return elements[element]; }
		const R& operator [] (int element) const { return elements[element]; }

		Vektor<R> getSubVektor(int startElement, int length);
		Vektor<R> getSubVektor(const std::vector<int>& retainedDimensions);

		R getGCD() const;
		Vektor<R> getGCDNormalizedVektor() const;

		static Vektor<R> getUnitVektor(int dim, int unitDim);
		static Vektor<R> getConstantVektor(int dim, R constant);
		static Vektor<R> getZeroVektor(int dim) { return getConstantVektor(dim, 0); }

		string toString() const;
		string toStringHtml(bool horizontal = true) const;
		string toStringJavaScript() const;
		void print() const { printf((toString() + "\n").c_str()); };

		bool isZero() const;
	};
}

#include "../cute/cute.h"

namespace AlgoTrans {
	template <class R> bool Vektor<R>::isZero() const {
		bool result = true;

		for (int q = 0; q < getLength(); q++) result = result && (elements[q] == 0);

		return result;
	}

	template <class R> string Vektor<R>::toStringJavaScript() const {
		string r = "new Vector([";

		for (int q = 0; q < getLength(); q++) { if (q > 0) r += ",";
			r += elements[q].toString();
		}

		return r + "])";
	}

	template <class R> string Vektor<R>::toString() const {
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

	/*template <class R> string Vektor<R>::toStringLatex(bool horizontal) const {
		string result = "\\left[ \\begin{array}{";
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

		return result + "\\end{array} \\right] ";
	}*/

	template <class R> Vektor<R> operator + (const Vektor<R>& a, const Vektor<R>& b) {
		ASSERT(a.getLength() == b.getLength());

		Vektor<R> result = Vektor<R>();

		for (int q = 0; q < a.getLength(); q++) result.appendElement(a[q] + b[q]);

		return result;
	};

	template <class R> Vektor<R> operator * (const R& f, const Vektor<R>& v) {
		Vektor<R> result = Vektor<R>();

		for (int q = 0; q < v.getLength(); q++) result.appendElement(f * v[q]);

		return result;
	};

	template <class R> Vektor<R> operator - (const Vektor<R>& a, const Vektor<R>& b) {
		ASSERT(a.getLength() == b.getLength());

		Vektor<R> result = Vektor<R>(a.getLength());

		for (int q = 0; q < a.getLength(); q++) result[q] = (a[q] - b[q]);

		return result;
	};

	template <class R> Vektor<R> Vektor<R>::operator - () const {
		Vektor<R> result = Vektor<R>(getLength());

		for (int q = 0; q < getLength(); q++) result[q] = -elements[q];

		return result;
	};

	template <class R> Vektor<R> operator << (const Vektor<R>& a, const Vektor<R>& b) {
		Vektor<R> result = Vektor<R>();

		for (int q = 0; q < a.getLength(); q++) result.appendElement(a[q]);
		for (int q = 0; q < b.getLength(); q++) result.appendElement(b[q]);

		return result;
	};

	template <class R> bool operator == (const Vektor<R>& a, const Vektor<R>& b) {
		bool equal = (a.getLength() == b.getLength());

		for (int e = 0; equal && (e < a.getLength()); e++) {
			equal = equal && (a[e] == b[e]);
		}

		return equal;
	};

	template <class R> bool operator != (const Vektor<R>& a, const Vektor<R>& b) { return !(a == b); }

	template <class R> Vektor<R> Vektor<R>::getUnitVektor(int dim, int unitDim) {
		Vektor<R> result = Vektor<R>();

		for (int q = 0; q < dim; q++) result.elements.push_back((R) ((q == unitDim) ? 1 : 0));

		return result;
	}

	template <class R> Vektor<R> Vektor<R>::getConstantVektor(int dim, R constant) {
		Vektor<R> result = Vektor<R>(dim);

		for (int q = 0; q < dim; q++) result[q] = constant;

		return result;
	}

	template <class R> Vektor<R>& Vektor<R>::operator << (const R& element) {
		appendElement(element);

		return *this;
	}

	template <class G> Vektor<G> operator << (const G& a, const Vektor<G>& b) { // hor concatenation
		Vektor<G> result;

		result.appendElement(a);
		for (int q = 0; q < b.getLength(); q++) result.appendElement(b[q]);

		return result;
	}

	template <class G> bool operator < (const Vektor<G>& a, const Vektor<G>& b) { // lexicographical ordering
		for (int q = 0; q < b.getLength(); q++) {
			if (a[q] < b[q]) { return true; }
			else if (a[q] > b[q]) { return false; }
		}

		return false;
	}

	template <class G> bool operator <= (const Vektor<G>& a, const Vektor<G>& b) { // lexicographical ordering
		for (int q = 0; q < b.getLength(); q++) {
			if (a[q] < b[q]) { return true; }
			else if (a[q] > b[q]) { return false; }
		}

		return true;
	}

	template <class R> Vektor<R> Vektor<R>::getSubVektor(int startElement, int length) {
		ASSERT((0 <= startElement) && (startElement + length <= getLength()));

		Vektor<R> result = Vektor<R>();

		for (int q = startElement; q < startElement + length; q++) result.elements.push_back((*this)[q]);

		return result;
	}

	template <class R> Vektor<R> Vektor<R>::getSubVektor(const std::vector<int>& retainedDimensions) {
		typedef std::vector<int>::const_iterator DimIterator;

		Vektor<R> result = Vektor<R>();

		for (DimIterator q = retainedDimensions.begin(); q != retainedDimensions.end(); ++q) {
			ASSERT(((*q) >= 0) && ((unsigned int)  (*q) < elements.size()));
			result.elements.push_back((*this)[*q]);
		}

		return result;
	}

	template <class R> Vektor<R>& Vektor<R>::operator += (const Vektor<R>& other) {
		for (int q = 0; q < getLength(); q++) elements[q] += other[q];

		return *this;
	}

	template <class R>
	std::string Vektor<R>::toStringHtml(bool horizontal) const {
		string result = "<table><tr><td>V</td><td>";
		if ((elements.size() != 0)) {
			result += "<table cellspacing=\"0\" border=\"1\">";
			if (horizontal) result += "<tr>";
			for (unsigned int r = 0; r < elements.size(); r++) {
				if (!horizontal) result += "<tr>";
				result += "<td align=\"center\" width=\"16\">" + elements[r].toString() + "</td>";
				if (!horizontal) result += "</tr>";
			}
			if (horizontal) result += "</tr>";
			result += "</table>";
		}

		return result + "</td></tr></table>";
	}

	template <class R>
	R Vektor<R>::getGCD() const {
		R g = 0;

		for (int q = 0; q < getLength(); q++) {
			g = gcd(elements[q], g);
		}

		return g;
	}

	template <class R>
	Vektor<R> Vektor<R>::getGCDNormalizedVektor() const {
		Vektor<R> result = Vektor<R>(getLength());

		R gcd = getGCD();
		if (gcd != 0) {
			for (int q = 0; q < getLength(); q++) result[q] = elements[q] / gcd;
		}

		return result;
	}
}

#endif /*CVECT_H_*/
