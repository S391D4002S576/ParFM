#ifndef CSPARSEVEKTOR_H_
#define CSPARSEVEKTOR_H_

#include <vector>
using std::vector;

#include <string>
using std::string;

#include "scalar/Bool.h"

namespace AlgoTrans {
	template <class R> class CSparseVector {
		template <class G> friend CSparseVector<G> operator << (const CSparseVector<G>& a, const CSparseVector<G>& b);
		template <class G> friend bool operator == (const CSparseVector<G>& a, const CSparseVector<G>& b);
	public:
	    struct DataPair {
	    	int dim;
	    	R value;

	    	DataPair(int iDim, R iValue) : dim(iDim), value(iValue) { }
	    };

	protected:
		R zeroConst;

		int dimension;

	    bool compareDataPair(const DataPair& a, const DataPair& b) { return a.dim < b.dim; }

		vector<DataPair> dataPairs;

		void sortDataPairs() { if (!sorted) sort(dataPairs.begin(), dataPairs.end(), compareDataPair); sorted = true; };
	public:
		bool sorted;

		void appendDataPair(int iDim, R iValue) { sorted = sorted && ((dataPairs.size() == 0) || (dataPairs[dataPairs.size() - 1].dim < iDim)); dataPairs.push_back(DataPair(iDim, iValue)); }
		void appendDataPair(DataPair iDataPair) { dataPairs.push_back(iDataPair); sorted = false; }
	public:
		vector<DataPair>& getDataPairs() { return dataPairs; }

		CSparseVector() : zeroConst(0), dimension(0), sorted(true) { };
		CSparseVector(int iDimension, vector<DataPair> iDataPairs) : zeroConst(0), dimension(iDimension), dataPairs(iDataPairs) {  }

		void appendElement(R element) { if (element != (R) 0) dataPairs.push_back(DataPair(dimension, element)); dimension++; }

		CSparseVector(int count, R iElements[]) : zeroConst(0) { dataPairs = vector<DataPair>(); for (int i = 0; i < count; i++) appendElement(iElements[i]); };

		CSparseVector(int count) : zeroConst(0), dimension(count) { dataPairs = vector<DataPair>(); };
		CSparseVector(int count, R defaultValue) : zeroConst(0), dimension(0) { dataPairs = vector<DataPair>(); if (defaultValue != (R) 0) for (int i = 0; i < count; i++) { appendElement(defaultValue); } else { dimension = count; } };

		CSparseVector(const CSparseVector<R>& iOriginal) : zeroConst(0) { if (this != &iOriginal) *this = iOriginal; };

		CSparseVector<R>& operator = (const CSparseVector<R>& other) {
			if (this != &other) {
				zeroConst = other.zeroConst;
				dimension = other.dimension;
				dataPairs = other.dataPairs;
				sorted = other.sorted;
			}

			return *this;
		}

		CSparseVector<R>& operator += (const CSparseVector<R>& other);

		CSparseVector<R>& operator << (const R& element) { appendElement(element); return *this; }

		CSparseVector<R> operator -();

		int getLength() const { return dimension; };

		R& operator [] (int element);
		const R& operator [] (int element) const;

		CSparseVector<R> getSubVector(int startElement, int length);
		CSparseVector<R> getSubVector(std::vector<int> retainedDimensions);

		R getGGD() const;

		static CSparseVector<R> getUnitVector(int dim, int unitDim);
		static CSparseVector<R> getConstantVector(int dim, R constant) { return CSparseVector<R>(dim, constant); }
		static CSparseVector<R> getZeroVector(int dim) { return CSparseVector(dim); }

		string toString(bool sparse = false) const;
		string toStringHtml(bool horizontal = true) const;
		void print() const { printf((toString() + "\n").c_str()); };

		bool isZero() const { return dataPairs.size() == 0; }
	};
}

#include "../cute/cute.h"

namespace AlgoTrans {
	template <class R> string CSparseVector<R>::toString(bool sparse) const {
		if (!sparse) {
			int maxL = 0;
			string strings[getLength()];
			for (int q = 0; q < getLength(); q++) {
				strings[q] = (*this)[q].toString();
				int strL = strings[q].length();
				maxL = (strL > maxL) ? strL : maxL;
			}

			string result;
			for (int q = 0; q < getLength(); q++) {
				for (int s = (maxL - strings[q].length()); s >= 0; s--) result += " ";
				result += strings[q];
			}

			return result;
		} else {
			string result;
			for (unsigned int d = 0; d < dataPairs.size(); d++) {
				if (d > 0) result += " ";
				result += CInteger(dataPairs[d].dim).toString() + ":" + dataPairs[d].value.toString();
			}
			return "Sparse(" + result + ")";
		}
	}

	/*template <class R> string CSparseVector<R>::toStringLatex(bool horizontal) const {
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

	template <class R> CSparseVector<R> operator + (const CSparseVector<R>& a, const CSparseVector<R>& b) {
		ASSERT(a.getLength() == b.getLength());
		ASSERT(a.sorted && b.sorted);

		CSparseVector<R> result = CSparseVector<R>(a.getLength());

		// Merge add
		int i = 0, j = 0;
		if (i == a.dataPairs.size()) {
			for (int q = j; q < b.dataPairs.size(); q++) result.addDataPair(b.dataPairs[q]);
		} else if (j == b.dataPairs.size()) {
			for (int q = i; q < a.dataPairs.size(); q++) result.addDataPair(a.dataPairs[q]);
		} else {
			if (a.dataPairs[i].dim == b.dataPairs[j].dim) {
				R sum = a.dataPairs[i] + b.dataPairs[j];
				if (sum != 0) result.addDataPair(a.dataPairs[i].dim, sum);
				i++; j++;
			} else if (a.dataPairs[i].dim < b.dataPairs[j].dim) { result.addDataPair(a.dataPairs[i++]);
			} else if (a.dataPairs[i].dim > b.dataPairs[j].dim) { result.addDataPair(b.dataPairs[j++]); }
		}

		result.sorted = true;

		return result;
	};

	template <class R> CSparseVector<R>& CSparseVector<R>::operator += (const CSparseVector<R>& other) {
		ASSERT(getLength() == other.getLength());
		ASSERT(sorted && other.sorted);

		CSparseVector<R> result = CSparseVector<R>(getLength());

		// Merge add
		unsigned int i = 0, j = 0;
		if (i == dataPairs.size()) {
			for (unsigned int q = j; q < other.dataPairs.size(); q++) result.appendDataPair(other.dataPairs[q]);
		} else if (j == other.dataPairs.size()) {
			for (unsigned int q = i; q < dataPairs.size(); q++) result.appendDataPair(dataPairs[q]);
		} else {
			if (dataPairs[i].dim == other.dataPairs[j].dim) {
				R sum = dataPairs[i].value + other.dataPairs[j].value;
				if (sum != (R) 0) result.appendDataPair(dataPairs[i].dim, sum);
				i++; j++;
			} else if (dataPairs[i].dim < other.dataPairs[j].dim) { result.appendDataPair(dataPairs[i++]);
			} else if (dataPairs[i].dim > other.dataPairs[j].dim) { result.appendDataPair(other.dataPairs[j++]); }
		}

		result.sorted = true;
		*this = result;

		return *this;
	};

	template <class R> CSparseVector<R> operator * (const R& f, const CSparseVector<R>& v) {
		ASSERT(v.sorted);

		vector<R> result = vector<R>();

		for (int q = 0; q < v.dataPairs.size(); q++) result.appendDataPair(v.dataPairs[q].dim, f * v.dataPairs[q].value);

		result.sorted = true;

		return result;
	};

	template <class R> CSparseVector<R> operator - (const CSparseVector<R>& a, const CSparseVector<R>& b) {
		ASSERT(a.getLength() == b.getLength());
		ASSERT(a.sorted && b.sorted);

		CSparseVector<R> result = CSparseVector<R>(a.getLength());

		// Merge add
		int i = 0, j = 0;
		if (i == a.dataPairs.size()) {
			for (int q = j; q < b.dataPairs.size(); q++) result.addDataPair(b.dataPairs[q]);
		} else if (j == b.dataPairs.size()) {
			for (int q = i; q < a.dataPairs.size(); q++) result.addDataPair(a.dataPairs[q]);
		} else {
			if (a.dataPairs[i].dim == b.dataPairs[j].dim) {
				R diff = a.dataPairs[i] - b.dataPairs[j];
				if (diff != 0) result.addDataPair(a.dataPairs[i].dim, diff);
				i++; j++;
			} else if (a.dataPairs[i].dim < b.dataPairs[j].dim) { result.addDataPair(a.dataPairs[i++]);
			} else if (a.dataPairs[i].dim > b.dataPairs[j].dim) { result.addDataPair(b.dataPairs[j++]); }
		}

		result.sorted = true;

		return result;
	};

	template <class R> CSparseVector<R> CSparseVector<R>::operator - () {
		ASSERT(sorted);

		CSparseVector<R> result = CSparseVector<R>(getLength());

		for (int q = 0; q < dataPairs.size(); q++) result.appendDataPair(dataPairs[q].dim, -dataPairs[q].value);

		result.sorted = true;

		return result;
	};

	template <class R> CSparseVector<R> operator << (const CSparseVector<R>& a, const CSparseVector<R>& b) {
		ASSERT(a.sorted);
		ASSERT(b.sorted);

		CSparseVector<R> result = CSparseVector<R>();

		result = a;
		for (unsigned int q = 0; q < b.dataPairs.size(); q++) result.appendDataPair(a.getLength() + b.dataPairs[q].dim, b.dataPairs[q].value);

		result.sorted = true;

		return result;
	};

	template <class R> bool operator == (const CSparseVector<R>& a, const CSparseVector<R>& b) {
		bool equal = (a.getLength() == b.getLength());

		for (unsigned int e = 0; equal && (e < a.dataPairs.size()); e++) {
			if (e >= b.dataPairs.size()) return false;

			equal = equal && (a[e] == b[e]);
		}

		return equal;
	};

	template <class R> bool operator != (const CSparseVector<R>& a, const CSparseVector<R>& b) { return !(a == b); }

	template <class R> CSparseVector<R> CSparseVector<R>::getUnitVector(int dim, int unitDim) {
		CSparseVector<R> result = CSparseVector<R>(dim);

		result.appendDataPair(unitDim, 1);

		result.sorted = true;

		return result;
	}

	template <class G> CSparseVector<G> operator << (const G& a, const CSparseVector<G>& b) { // hor concatenation
		CSparseVector<G> result;

		result.appendElement(a);

		return result << b;
	}

	template <class G> bool operator < (const CSparseVector<G>& a, const CSparseVector<G>& b) { // lexicographical ordering
		// Merge compare
		int i = 0, j = 0;
		if (i == a.dataPairs.size()) {
			for (int q = j; q < b.dataPairs.size(); q++) if (b.dataPairs[j].value > 0) return true;
		} else if (j == b.dataPairs.size()) {
			for (int q = i; q < a.dataPairs.size(); q++) if (a.dataPairs[i].value > 0) return false;
		} else {
			if (a.dataPairs[i].dim == b.dataPairs[j].dim) {
				if (a.dataPairs[i].value < b.dataPairs[j].value) return false;
				if (a.dataPairs[i].value > b.dataPairs[j].value) return true;
				i++; j++;
			} else if (a.dataPairs[i].dim < b.dataPairs[j].dim) { if (a.dataPairs[i].value > 0) return false;
			} else if (a.dataPairs[i].dim > b.dataPairs[j].dim) { if (b.dataPairs[j].value > 0) return true; }
		}

		return false;
	}

	template <class R> CSparseVector<R> CSparseVector<R>::getSubVector(int startElement, int length) {
		ASSERT((0 <= startElement) && (startElement + length <= getLength()));

		CSparseVector<R> result = CSparseVector<R>(length);

		for (int q = 0; q < dataPairs.size(); q++) {
			if ((dataPairs[q].dim >= startElement) && (dataPairs[q].dim < startElement + length)) {
				result.appendDataPair(dataPairs[q].dim - startElement, dataPairs[q].value);
			}
		}

		return result;
	}

	template <class R> CSparseVector<R> CSparseVector<R>::getSubVector(std::vector<int> retainedDimensions) {
		typedef std::vector<int>::const_iterator DimIterator;

		CSparseVector<R> result = CSparseVector<R>();

		std::sort(retainedDimensions.begin(), retainedDimensions.end());
		unsigned int rd = 0;
		for (unsigned int q = 0; q < dataPairs.size(); q++) {
			while ((rd < retainedDimensions.size()) && (retainedDimensions[rd] < dataPairs[q].dim)) rd++;
			if (rd >= retainedDimensions.size()) break;

			if (dataPairs[q].dim == retainedDimensions[rd]) result.appendDataPair(rd, dataPairs[q].value);
		}

		result.sorted = true;

		return result;
	}

	template <class R> R CSparseVector<R>::getGGD() const {
		R result = 1;

		for (int q = 0; q < dataPairs.size(); q++) result = ggd(result, dataPairs[q].value);

		return (dataPairs.size() > 0) ? result : 0;
	}

	template <class R> R& CSparseVector<R>::operator [] (int element) {
		// XXX: replace with binary search
		for (unsigned int q = 0; q < dataPairs.size(); q++) {
			if (dataPairs[q].dim == element) return dataPairs[q].value;
			if (dataPairs[q].dim > element) break;
		}

		// Replace this vector with an extended vector with a 0 in the right position... (messy!)
		bool inserted = false;
		CSparseVector<R> resultV = CSparseVector<R>(getLength());
		for (unsigned int q = 0; q < dataPairs.size(); q++) {
			if (!inserted && (dataPairs[q].dim > element)) {
				resultV.appendDataPair(element, 0);
				inserted = true;
			}
			resultV.appendDataPair(dataPairs[q]);
		}
		if (!inserted) resultV.appendDataPair(element, 0);

		*this = resultV;
		sorted = true;

		return (*this)[element];
	}

	template <class R> const R& CSparseVector<R>::operator [] (int element) const {
		// XXX: replace with binary search
		for (unsigned int q = 0; q < dataPairs.size(); q++) {
			if (dataPairs[q].dim == element) return dataPairs[q].value;
			if (dataPairs[q].dim > element) return zeroConst;
		}

		return zeroConst;
	}

	template <class R>
	std::string CSparseVector<R>::toStringHtml(bool horizontal) const {
		string result = "<table><tr><td>Vec</td><td>";
		if ((getLength() != 0)) {
			result += "<table cellspacing=\"0\" border=\"1\">";
			if (horizontal) result += "<tr>";
			for (unsigned int r = 0; r < getLength(); r++) {
				if (!horizontal) result += "<tr>";
				result += "<td>" + (*this)[r].toString() + "</td>";
				if (!horizontal) result += "</tr>";
			}
			if (horizontal) result += "</tr>";
			result += "</table>";
		}

		return result + "</td></tr></table>";
	}
}

#endif /*CVECT_H_*/
