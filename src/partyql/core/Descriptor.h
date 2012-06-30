#ifndef DESKRIPTOR_H_
#define DESKRIPTOR_H_

#include <vector>
#include <list>
#include <string>
#include <algorithm>

#include "Config.h"

#include "Declarations.h"

#include "scalar/Integer.h"
#include "Vektor.h"
//#include "SparseVektor.h"
#include "util/Html.h"
#include "util/Util.h"

#include "IntegerProgramming.h"

#include "HadronDescription.h"

#include <glpk.h>

namespace AlgoTrans {
	template <class D> class DeskriptorSet;

	template <class R> class Deskriptor {
		template <class G> friend bool operator < (const Deskriptor<G>& a, const Deskriptor<G>& b);
		template <class G> friend bool operator == (const Deskriptor<G>& a, const Deskriptor<G>& b);
		template <class G> friend bool operator != (const Deskriptor<G>& a, const Deskriptor<G>& b);

		// Vector based:
		template <class G> friend G operator % (const Deskriptor<G>& a, const Deskriptor<G>& b);
		template <class G> friend Deskriptor<G> operator * (G& a, const Deskriptor<G>& b);
		template <class G> friend Deskriptor<G> operator + (const Deskriptor<G>& a, const Deskriptor<G>& b);
		template <class G> friend Deskriptor<G> operator - (const Deskriptor<G>& a, const Deskriptor<G>& b);
		template <class G> friend Deskriptor<G> operator << (Deskriptor<G>& a, Deskriptor<G>& b);
	public:
		typedef Vektor<R> Vector;
		typedef R DeskriptorElementType;
	private:
		Vector vector;
		bool bidir;
	public:
		Deskriptor() : vector(Vector()), bidir(false) { }
		Deskriptor(const Vector& iVector, bool iBidir) : vector(iVector), bidir(iBidir) { }

		bool isBidirectional() const { return bidir; }

		// Vector-based:
		Vector& getVector() { return vector; }
		const Vector& getVector() const { return vector; }

		R& operator [] (int i) { return vector[i]; };
		const R& operator [] (int i) const { return vector[i]; };
		const R& operator () (int i) const { return vector[i]; }; // forced const (why is this necessary?)

		int getSpaceDimension() const { return vector.getLength(); }

		R getGCD() const;
		Deskriptor<R> getDimensionPermutation(std::vector<int> sourceDimensions) const;
		Deskriptor<R> getSimplifiedDeskriptor() const;

		Deskriptor<R> operator -() const { return Deskriptor<R>(-vector, isBidirectional()); }; // invers op
		Deskriptor& operator += (const Deskriptor& rhs) { vector += rhs.getVector(); return *this; };
		Deskriptor<R> getProjection(std::vector<int> retainedDimensions) { return Deskriptor<R>(vector.getSubVektor(retainedDimensions), isBidirectional()); }
		Deskriptor<R> getProjection(int startDim, int dimCount) { return getProjection(getInterval(startDim, dimCount)); }
		Deskriptor<R> getProjection(const DeskriptorSet<Deskriptor<R> >& orthSpaceGenerators, bool preOrthogonalize = true) const;

		static Deskriptor<R> getUnitDeskriptor(int dimension, int unitDimension, bool iBidir) { return Deskriptor<R>(Vector::getUnitVektor(dimension, unitDimension), iBidir); };
		static Deskriptor<R> getConstantDeskriptor(int dimension, R constant, bool iBidir) { return Deskriptor<R>(Vector::getConstantVektor(dimension, constant), iBidir); };
		static Deskriptor<R> getZeroDeskriptor(int dimension, bool iBidir) { return Deskriptor<R>(Vector::getZeroVektor(dimension), iBidir); };

		// Erroneous?:
		//static Deskriptor<R> getTransposedDeskriptor(const DeskriptorSet<Deskriptor<R> >& deskriptors, int dimension);

		// Text
		std::string toString() const { return getTypeName() + "(" + vector.toString() + (isBidirectional() ? " == " : ">=") + ")"; }
		std::string toStringHtml() const { return Html::wrapShort(getTypeNameHtml(), "<table><tr><td>" + vector.toStringHtml() + "</td><td>" + (isBidirectional() ? " == " : ">=") + "</td></tr></table>"); }
		std::string toStringJavaScript() const;
		void print() const { printf("%s\n", toString().c_str()); };

		static std::string getTypeName() { return (SHORT_TYPENAMES ? "D<" : "Deskriptor<" ) + R::getTypeName() + ">"; }
		static std::string getTypeNameHtml() { return Html::wrapType((SHORT_TYPENAMES ? "D" : "Deskriptor" ), R::getTypeNameHtml()); }
	};

	template <class G> Deskriptor<G> operator + (const Deskriptor<G>& a, const Deskriptor<G>& b) { ASSERT(a.isBidirectional() == b.isBidirectional());
		return Deskriptor<G>(a.getVector() + b.getVector(), a.isBidirectional());
	}

	template <class G> Deskriptor<G> operator - (const Deskriptor<G>& a, const Deskriptor<G>& b) { 	return a + (-b); }

	template <class G> Deskriptor<G> operator * (G& a, const Deskriptor<G>& b) {
		return Deskriptor<G>(a * b.getVector(), b.isBidirectional());
	}

	template <class G> Deskriptor<G> operator << (Deskriptor<G>& a, Deskriptor<G>& b) { ASSERT(a.isBidirectional() == b.isBidirectional());
		return Deskriptor<G>(a.getVector() << b.getVector(), a.isBidirectional());
	}

	template <class G>
	std::string Deskriptor<G>::toStringJavaScript() const {
		return std::string("new Deskriptor(") + (isBidirectional() ? "true" : "false") + ", " + vector.toStringJavaScript() + ")";
	}

	template <class G> bool operator < (const Deskriptor<G>& a, const Deskriptor<G>& b) {
		if (a.isBidirectional() != b.isBidirectional()) return a.isBidirectional();

		return (a.getVector() < b.getVector());
	}

	template <class G> bool operator == (const Deskriptor<G>& a, const Deskriptor<G>& b) { return (a.bidir == b.bidir) && (a.vector == b.vector); }
	template <class G> bool operator != (const Deskriptor<G>& a, const Deskriptor<G>& b) { return !(a == b); }
	//template <class R> bool pVectorCompare(Deskriptor<R>* d1, Deskriptor<R>* d2) { return *d2 < *d1; }

	template <class R> bool pDeskriptorCompare(const Deskriptor<R>* d1, const Deskriptor<R>* d2) { return (*d1) < (*d2); }

	/*
	 * Deskriptors must have:
	 * - Equality test
	 *   - Vector: Equality per element
	 * - Total sortability --> comparision op
	 *   - Vector: Lexicographic ord
	 */
	template <class D> class DeskriptorSet {
		// Vector based:
		template <class DD> friend DeskriptorSet<DD> operator << (DeskriptorSet<DD>& a, DD& b);
		template <class DD> friend DeskriptorSet<DD> operator << (DeskriptorSet<DD> a, DeskriptorSet<DD> b);
		template <class DD> friend DeskriptorSet<DD> operator >> (const DeskriptorSet<DD>& a, const DeskriptorSet<DD>& b);
		template <class DD> friend DeskriptorSet<DD> operator >> (const DeskriptorSet<DD>& a, const DD& b);
	private:
		typedef typename D::DeskriptorElementType R;

		int dimension;
		void cleanup();

		bool minimal; // to avoid minimizing when unnecessary
		void actualMinimisation();

		bool sorted;
		void actualSort();

		bool bidNormalized;
		int bidirDeskCount;
		void normalizeBidirectionalDeskriptors();
		void actualNormalizeBidirectionalDeskriptors();

		DeskriptorSet getUnidirSubSet() const;

		// for bidir minimisation
		int rank;
		std::vector<int> pivotDims, nonPivotDims;
	protected:
		std::vector<D* > deskriptors;

		void swapDeskriptors(int ixA, int ixB) { D* d = deskriptors[ixA]; deskriptors[ixA] = deskriptors[ixB]; deskriptors[ixB] = d; }
	public:
		DeskriptorSet getBidirSubSet() const;
		void extractImplicitBidirectionalDeskriptors(bool retainOnlyBidirs = false);

		typedef D DeskriptorType;
		typedef DeskriptorSet* Pointer;
		//typedef vector<D*>::iterator iterator;
		//typedef vector<D*>::const_iterator const_iterator;

		//iterator begin() { return deskriptors.begin(); }
		//const_iterator begin() const { return deskriptors.begin(); }
		//iterator end() { return deskriptors.end(); }
		//const_iterator end() const { return deskriptors.end(); }

		//DeskriptorSet(int i) { data = (i == 1) ? Bool(true) : Bool(false); }

		DeskriptorSet(int iDimension) : dimension(iDimension), minimal(true), sorted(true), bidNormalized(true), bidirDeskCount(0) { };
		DeskriptorSet(const DeskriptorSet& iOriginal) { *this = iOriginal; }
		DeskriptorSet& operator = (const DeskriptorSet& other);

		~DeskriptorSet() { cleanup(); };

		D& operator [] (int i) { return *deskriptors[i]; minimal = false; sorted = false; };
		const D& operator [] (int i) const { return *deskriptors[i]; };

		int getSize() const { return deskriptors.size(); }
		int getDimension() const { return dimension; };
		int getBidirDeskriptorCount() const { return bidirDeskCount; }

		void addDeskriptor(D* deskriptor); // Ownership of deskriptor transferred to DeskriptorSet
		void addDeskriptor(const D& deskriptor) { addDeskriptor(new D(deskriptor)); minimal = false; sorted = false; bidNormalized = false; }; // Ownership of deskriptor remains with caller

		//DeskriptorSet getFMProjection(int projectionDim, bool preSimplify) const;

		DeskriptorSet getDual() const;

		DeskriptorSet getConstrainingProjection(std::vector<int> retainedDimensions) const;
		DeskriptorSet getFMProjection(int projectionDim, bool presimplify) const;
		DeskriptorSet getFMProjection_Deb(int projectionDim, bool presimplify) const;
		DeskriptorSet getGeneratingProjection(std::vector<int> retainedDimensions) const;
		DeskriptorSet getProjection(const DeskriptorSet& orthSpaceGenerators, bool preOrthogonalize = true) const;
		DeskriptorSet removeLinearGeneratedVertices(const DeskriptorSet& orthSpaceGenerators, bool preOrthogonalize = true) const;

		DeskriptorSet getOrthogonalSet() const;

		DeskriptorSet getSimplifiedDeskriptorSet(int numberOfCertainlyNonRedundant = 0) const;

		DeskriptorSet getTransposedDeskriptorSet(bool bidir) const;

		bool bidirectionallyGenerates(const D& deskriptor) const;
		bool nonNegativelyGenerates(const D& deskriptor) const;

		void minimize(); // XXX: Should find a way to prevent overriding this method, so actualMinimisation is overridden instead
		void sort(); // XXX: ...same

		static DeskriptorSet<D> unitDeskriptorSet(int dimension, bool bidir);
		static DeskriptorSet<D> zeroDeskriptorSet(int count, int dimension, bool bidir);
		static DeskriptorSet<D> fromSingleElement(D e);

		bool isEmpty() { return deskriptors.size() == 0; }

		// Vector-based
		DeskriptorSet<D> operator -() const;
		DeskriptorSet<D> getDimensionPermutation(std::vector<int> sourceDimensions) const;

		// Text
		std::string toString() const;
		std::string toStringHtml() const;
		std::string toStringJavaScript() const;

		void print(const std::string prefix = "") const { printf("%s: %s\n", prefix.c_str(), toString().c_str()); };

		static std::string getTypeName() { return (SHORT_TYPENAMES ? "DS<" : "DeskriptorSet<" ) + D::getTypeName() + ">"; }
		static std::string getTypeNameHtml() { return Html::wrapType((SHORT_TYPENAMES ? "DS" : "DeskriptorSet" ), D::getTypeNameHtml()); }
	};

	template <class D> DeskriptorSet<D> operator !(DeskriptorSet<D>& a) { return a.getDual(); }

	template <class R>
	Deskriptor<R> Deskriptor<R>::getDimensionPermutation(std::vector<int> sourceDimensions) const { ASSERT((int) sourceDimensions.size() == getSpaceDimension());

		Vector w = Vector();
		for (int q = 0; q < getSpaceDimension(); q++) w.appendElement(vector[sourceDimensions[q]]);

		return Deskriptor<R>(w, bidir);
	}

	template <class D>
	DeskriptorSet<D> DeskriptorSet<D>::fromSingleElement(D e) {
		DeskriptorSet<D> result = DeskriptorSet<D>(e.getSpaceDimension());

		result.addDeskriptor(e);

		return result;
	}

	template <class D>
	DeskriptorSet<D> DeskriptorSet<D>::getBidirSubSet() const {
		DeskriptorSet<D> result = DeskriptorSet<D>(getDimension());

		for (int q = 0; q < getSize(); q++) if ((*this)[q].isBidirectional()) result.addDeskriptor((*this)[q]);

		return result;
	}

	template <class D>
	DeskriptorSet<D> DeskriptorSet<D>::getOrthogonalSet() const {
		DeskriptorSet<D> result = DeskriptorSet<D>(getDimension());

		for (int q = 0; q < getSize(); q++) {
			result.addDeskriptor((*this)[q].getProjection(result, false));
		}
//		result.normalizeBidirectionalDeskriptors();
//		result.extractImplicitBidirectionalDeskriptors();

		return result;
	}

	template <class D>
	DeskriptorSet<D> DeskriptorSet<D>::getUnidirSubSet() const {
		DeskriptorSet<D> result = DeskriptorSet<D>(getDimension());

		for (int q = 0; q < getSize(); q++) if (!(*this)[q].isBidirectional()) result.addDeskriptor((*this)[q]);

		return result;
	}

	template <class D>
	DeskriptorSet<D> DeskriptorSet<D>::operator - () const {
		DeskriptorSet<D> result = DeskriptorSet<D>(getDimension());

		for (int q = 0; q < getSize(); q++) result.addDeskriptor(-((*this)[q]));

		return result;
	}

	template <class D>
	DeskriptorSet<D> DeskriptorSet<D>::getDimensionPermutation(std::vector<int> sourceDimensions) const {
		DeskriptorSet<D> result = DeskriptorSet<D>(getDimension());

		for (int q = 0; q < getSize(); q++) result.addDeskriptor((*this)[q].getDimensionPermutation(sourceDimensions));

		return result;
	}

	template <class R>
	Deskriptor<R> operator << (const Deskriptor<R>& a, const Deskriptor<R>& b) { ASSERT(a.isBidirectional() == b.isBidirectional()); // MMM
		return Deskriptor<R>(a.getVector() << b.getVector(), a.isBidirectional());
	}

	template <class DD> DeskriptorSet<DD> operator << (DeskriptorSet<DD>& a, DD& b) {
		DeskriptorSet<DD> result = DeskriptorSet<DD>();

		for (int q = 0; q < a.getDeskriptorCount(); q++) result.addDeskriptor(a[q] << b);

		return result;
	}

	template <class DD> DeskriptorSet<DD> operator << (DeskriptorSet<DD> a, DeskriptorSet<DD> b) {
		ASSERT(a.getSize() == b.getSize());

		DeskriptorSet<DD> result = DeskriptorSet<DD>(a.getDimension() + b.getDimension());

		for (int q = 0; q < a.getSize(); q++) result.addDeskriptor(a[q] << b[q]);

		return result;
	}

	//template <class DD> DeskriptorSet<DD> operator >> (const DeskriptorSet<DD>& a, const DeskriptorSet<DD>& b) { return a + b; }
	template <class DD> DeskriptorSet<DD> operator >> (const DeskriptorSet<DD>& a, const DD& b) { return a >> DeskriptorSet<DD>::fromSingleElement(b); }

	template <class D>
	DeskriptorSet<D> operator >> (const DeskriptorSet<D>& a, const DeskriptorSet<D>& b) {
		ASSERT(a.getDimension() == b.getDimension());

		DeskriptorSet<D> result = DeskriptorSet<D>(a.getDimension());

		for (int q = 0; q < a.getSize(); q++) result.addDeskriptor(a[q]);
		for (int q = 0; q < b.getSize(); q++) result.addDeskriptor(b[q]);

		//result.minimize(); // --> efficiency may depend on surrounding code, so this may not be ideal

		return result;
	}

	template <class D>
	DeskriptorSet<D> operator + (const DeskriptorSet<D>& a, const DeskriptorSet<D>& b) {
		ASSERT(a.getDimension() == b.getDimension());

		DeskriptorSet<D> result = DeskriptorSet<D>(a.getDimension());

		for (int q = 0; q < a.getSize(); q++) result.addDeskriptor(a[q]);
		for (int q = 0; q < b.getSize(); q++) result.addDeskriptor(b[q]);

		result.minimize(); // --> efficiency may depend on surrounding code, so this may not be ideal

		return result;
	}

	template <class D>
	DeskriptorSet<D> DeskriptorSet<D>::unitDeskriptorSet(int dimension, bool bidir) {
		DeskriptorSet<D> result = DeskriptorSet<D>(dimension);

		for (int q = 0; q < dimension; q++) result.addDeskriptor(D::getUnitDeskriptor(dimension, q, bidir));

		return result;
	}

	template <class R>
	R operator % (const Deskriptor<R>& a, const Deskriptor<R>& b) { ASSERT(a.getSpaceDimension() == b.getSpaceDimension());
		R sum = 0;

		for (int q = 0; q < a.getSpaceDimension(); q++) sum += a[q] * b[q];

		return sum;
	}

	template <class D>
	DeskriptorSet<D> DeskriptorSet<D>::zeroDeskriptorSet(int count, int dimension, bool bidir) {
		DeskriptorSet<D> result = DeskriptorSet<D>(dimension);

		for (int q = 0; q < count; q++) result.addDeskriptor(D::getZeroDeskriptor(dimension, bidir));

		return result;
	}

	template <class D>
	void DeskriptorSet<D>::cleanup() {
		for (unsigned int q = 0; q < deskriptors.size(); q++) delete deskriptors[q];
		deskriptors.clear();
	}

	template <class D>
	void DeskriptorSet<D>::minimize() {
		if (!minimal) {
			actualMinimisation();
			minimal = true;
		}
	}

	template <class D>
	void DeskriptorSet<D>::actualSort() {
		 std::sort(deskriptors.begin(), deskriptors.end(), pDeskriptorCompare<typename D::DeskriptorElementType>);
	}

	template <class D>
	void DeskriptorSet<D>::sort() {
		if (!sorted) {
			actualSort();
			sorted = true;
		}
	}

	/*template <class D>
	DeskriptorSet<D> DeskriptorSet<D>::getConstrainingProjection(std::vector<int> retainedDimensions) const {
		// XXX: Only works for bidir scubes (modifications required for other hadrons)
		// mods: cons with zero prefix may still constrain somewhat (???)
		std::vector<int> removedDimensions = getDimSetComplement(retainedDimensions, getDimension());
		printf("Projecting %d dimensions from %d dimensional set, retaining %d\n", (unsigned int) removedDimensions.size(), getDimension(), (unsigned int) retainedDimensions.size());

		// Reorder dims so minimization begins with removed dims
		DeskriptorSet<D> cons = DeskriptorSet<D>(getDimension());
		for (int q = 0; q < getSize(); q++) {
			DeskriptorType d = DeskriptorType(deskriptors[q]->getProjection(removedDimensions))
								<< DeskriptorType(deskriptors[q]->getProjection(retainedDimensions));
			cons.addDeskriptor(new DeskriptorType(d));
		}
		cons.minimize();

		DeskriptorSet<D> resultC = DeskriptorSet<D>(retainedDimensions.size());
		for (int q = 0; q < cons.getSize(); q++) {
			if (cons[q].getProjection(getInterval(0, removedDimensions.size())) == DeskriptorType::getZeroDeskriptor(removedDimensions.size())) {
				resultC.addDeskriptor(new DeskriptorType(cons[q].getProjection(getInterval(removedDimensions.size(), retainedDimensions.size()))));
			} // other constraints may still refine the domain described by the retained dimensions... XXX
		}

		return resultC;
	}*/

	template <class D>
	DeskriptorSet<D> DeskriptorSet<D>::getConstrainingProjection(std::vector<int> retainedDimensions) const {
		std::vector<int> removedDimensions = getDimSetComplement(retainedDimensions, getDimension());
		std::sort(removedDimensions.begin(), removedDimensions.end());

	    if (tcConstrainingProjection == NULL) { tcConstrainingProjection = new TimeCollector("Constraining Projection"); timeCollectors.push_back(tcConstrainingProjection); }
	    tcConstrainingProjection->resume();

		DeskriptorSet<D> result = *this;
		//result.print();
		for (int q = removedDimensions.size() - 1; q >= 0; q--) {
//			CInteger(removedDimensions[q]).print();
			result = result.getFMProjection(removedDimensions[q], q == ((int) removedDimensions.size()) - 1);
			//result.print("fm" + CInteger(q).toString());
			//result.print();
		}

	    tcConstrainingProjection->pause();

		return result;
	}

	template <class D>
	DeskriptorSet<D> DeskriptorSet<D>::getFMProjection(int projectionDim, bool preSimplify) const {
		typedef typename D::DeskriptorElementType R;
		ASSERT((projectionDim >= 0) && (projectionDim < getDimension()));

		DeskriptorSet<D> Y = preSimplify ? getSimplifiedDeskriptorSet() : *this;
		//Y.print();
		DeskriptorSet<D> result = DeskriptorSet<D>(getDimension() - 1);

		std::vector<D*> posList, negList; // Lists with constraints with zero-coeff, pos coeff, neg coeff for projectionDim respectively
		std::vector<D*> eqList;
		for (int r = 0; r < Y.getSize(); r++) {
			R pivot = Y[r][projectionDim];
			if (pivot == 0) {
  				  result.addDeskriptor(Y[r].getProjection(0, projectionDim)
									 << Y[r].getProjection(projectionDim + 1, result.getDimension() - projectionDim));
			} else {
				if (!Y[r].isBidirectional()) {
					if (pivot > 0) { posList.push_back(&Y[r]); }
					else { negList.push_back(&Y[r]); }
				} else {
					eqList.push_back(&Y[r]);
				}
			}
		}

		//int nonRedCount = result.getSize();

		// Recombine pos & neg
		//CInteger(posList.size()*negList.size());
		for (unsigned int p = 0; p < posList.size(); p++) { D& rP = *posList[p];
			for (unsigned int n = 0; n < negList.size(); n++) { D& rN = *negList[n];
				Vektor<R> row = Vektor<R>(result.getDimension());
				bool isZero = true;
				for (int c = 0; c < result.getDimension(); c++) { int d = (c < projectionDim) ? c : c + 1;
					row[c] = rP[projectionDim] * rN[d] - rP[d] * rN[projectionDim];
					isZero = isZero && (row[c] == 0);
				}
				if (!isZero) result.addDeskriptor(D(row, false));
			}
		}

		for (unsigned int p = 0; p < posList.size(); p++) { D& rP = *posList[p];
			for (unsigned int n = 0; n < eqList.size(); n++) { D& rE = *eqList[n];
				Vektor<R> row = Vektor<R>(result.getDimension());
			    bool isZero = true;
				if (rE[projectionDim] < 0) {
				    for (int c = 0; c < result.getDimension(); c++) { int d = (c < projectionDim) ? c : c + 1;
  					    row[c] = rP[projectionDim] * rE[d] - rP[d] * rE[projectionDim];
					    isZero = isZero && (row[c] == 0);
				    }
				} else if (rE[projectionDim] > 0) {
				    for (int c = 0; c < result.getDimension(); c++) { int d = (c < projectionDim) ? c : c + 1;
  					    row[c] = -rP[projectionDim] * rE[d] + rP[d] * rE[projectionDim];
					    isZero = isZero && (row[c] == 0);
				    }
				}
				if (!isZero) result.addDeskriptor(D(row, false));
			}
		}

		for (unsigned int p = 0; p < eqList.size(); p++) { D& rE = *eqList[p];
			for (unsigned int n = 0; n < negList.size(); n++) { D& rN = *negList[n];
				Vektor<R> row = Vektor<R>(result.getDimension());
				    bool isZero = true;
				if (rE[projectionDim] > 0) {
				    for (int c = 0; c < result.getDimension(); c++) { int d = (c < projectionDim) ? c : c + 1;
  					    row[c] = rE[projectionDim] * rN[d] - rE[d] * rN[projectionDim];
					    isZero = isZero && (row[c] == 0);
				    }
				} else {
				    for (int c = 0; c < result.getDimension(); c++) { int d = (c < projectionDim) ? c : c + 1;
  					    row[c] = -rE[projectionDim] * rN[d] + rE[d] * rN[projectionDim];
					    isZero = isZero && (row[c] == 0);
				    }
				}
				if (!isZero) result.addDeskriptor(D(row, false));
			}
		}

		for (unsigned int p = 0; p < eqList.size(); p++) { D& rE = *eqList[p];
			for (unsigned int n = 0; n < eqList.size(); n++) { D& rN = *eqList[n];
				Vektor<R> row = Vektor<R>(result.getDimension());
				bool isZero = true;
			    //if ((rE[projectionDim] > 0) && (rN[projectionDim] < 0)) {
				    for (int c = 0; c < result.getDimension(); c++) { int d = (c < projectionDim) ? c : c + 1;
  					    row[c] = rE[projectionDim] * rN[d] - rE[d] * rN[projectionDim];
  					    isZero = isZero && (row[c] == 0);
				    }
		        /*} else if ((rN[projectionDim] > 0) && (rE[projectionDim] < 0)) {
				    for (int c = 0; c < result.getDimension(); c++) { int d = (c < projectionDim) ? c : c + 1;
  					    row[c] = rN[projectionDim] * rE[d] - rN[d] * rE[projectionDim];
  					    isZero = isZero && (row[c] == 0);
				    }
		        } else if ((rN[projectionDim] > 0) && (rE[projectionDim] > 0)) {
				    for (int c = 0; c < result.getDimension(); c++) { int d = (c < projectionDim) ? c : c + 1;
  					    row[c] = rN[projectionDim] * rE[d] + rN[d] * rE[projectionDim];
  					    isZero = isZero && (row[c] == 0);
				    }
		        } else if ((rN[projectionDim] < 0) && (rE[projectionDim] < 0)) {
				    for (int c = 0; c < result.getDimension(); c++) { int d = (c < projectionDim) ? c : c + 1;
  					    row[c] = -rN[projectionDim] * rE[d] - rN[d] * rE[projectionDim];
  					    isZero = isZero && (row[c] == 0);
				    }
		        }*/
			    if (!isZero) result.addDeskriptor(D(row, true));
			}
		} // XXX: inefficient

		ASSERT(result.bidirDeskCount >= 0);

		return result.getSimplifiedDeskriptorSet(0); //nonRedCount);
	}

	template <class D>
	DeskriptorSet<D> DeskriptorSet<D>::getFMProjection_Deb(int projectionDim, bool preSimplify) const {
		typedef typename D::DeskriptorElementType R;
		ASSERT((projectionDim >= 0) && (projectionDim < getDimension()));

		DeskriptorSet<D> Y = preSimplify ? getSimplifiedDeskriptorSet() : *this;
		DeskriptorSet<D> Z = DeskriptorSet<D>(getDimension());
		//Y.print();
		DeskriptorSet<D> result = DeskriptorSet<D>(getDimension() - 1);

		std::vector<D*> posList, negList; // Lists with constraints with zero-coeff, pos coeff, neg coeff for projectionDim respectively
		std::vector<D*> eqList;
		for (int r = 0; r < Y.getSize(); r++) {
			R pivot = Y[r][projectionDim];
			if (pivot == 0) {
				result.addDeskriptor(Y[r].getProjection(0, projectionDim)
									 << Y[r].getProjection(projectionDim + 1, result.getDimension() - projectionDim));
			} else {
				if (!Y[r].isBidirectional()) {
					if (pivot > 0) { posList.push_back(&Y[r]); }
					else { negList.push_back(&Y[r]); }
				} else {
					eqList.push_back(&Y[r]);
					if (pivot > 0) {
					  posList.push_back(&Y[r]);
					  Z.addDeskriptor(-Y[r]);
					  negList.push_back(&Z[Z.getSize() - 1]);
					} else {
						negList.push_back(&Y[r]);
						Z.addDeskriptor(-Y[r]);
						posList.push_back(&Z[Z.getSize() - 1]);
					}
				}
			}
		}

		int nonRedCount = result.getSize();

		// Recombine pos & neg
		//CInteger(posList.size()*negList.size());
		for (unsigned int p = 0; p < posList.size(); p++) { D& rP = *posList[p];
			for (unsigned int n = 0; n < negList.size(); n++) { D& rN = *negList[n];
				Vektor<R> row = Vektor<R>(result.getDimension());
				bool isZero = true;
				for (int c = 0; c < result.getDimension(); c++) { int d = (c < projectionDim) ? c : c + 1;
					row[c] = rP[projectionDim] * rN[d] - rP[d] * rN[projectionDim];
					isZero = isZero && (row[c] == 0);
				}
				if (!isZero) result.addDeskriptor(D(row, false));
			}
		}

		//result.print();
		return result.getSimplifiedDeskriptorSet(nonRedCount);
	}

	template <class D>
	DeskriptorSet<D> DeskriptorSet<D>::getGeneratingProjection(std::vector<int> retainedDimensions) const {
		// XXX: Only works for vectorized sets (modifications required for other hadrons)
		DeskriptorSet<D> resultG = DeskriptorSet<D>(retainedDimensions.size());
		for (int q = 0; q < getSize(); q++) {
			resultG.addDeskriptor(new DeskriptorType(deskriptors[q]->getProjection(retainedDimensions)));
		}

		resultG.minimize();
		ASSERT(resultG.bidirDeskCount >= 0);

		return resultG;
	}

	template <class D> //--> almost original version for non-sparse vektors
	DeskriptorSet<D> DeskriptorSet<D>::getDual() const {
		DeskriptorSet<D> d = -unitDeskriptorSet(getDimension(), true);
		DeskriptorSet<D> combSet = DeskriptorSet<D>(getDimension() + getSize());
		DeskriptorSet<D> transpSet = getTransposedDeskriptorSet(true);

		if (tcDualisation == NULL) { tcDualisation = new TimeCollector("Dualisation"); timeCollectors.push_back(tcDualisation); }
		tcDualisation->resume();

		//this->print();
		std::vector<int> retainedDimensions = std::vector<int>();
		for (int q = 0; q < getDimension(); q++) {
			combSet.addDeskriptor(d[q] << transpSet[q]);
			retainedDimensions.push_back(q);
		}
		for (int q = 0; q < getSize(); q++) {
		    if (!deskriptors[q]->isBidirectional()) combSet.addDeskriptor(D::getZeroDeskriptor(getDimension(), false) << D::getUnitDeskriptor(getSize(), q, false));
		}
		//combSet = combSet >> (DeskriptorSet<D>::zeroDeskriptorSet(getSize(), getDimension(), false) << unitDeskriptorSet(getSize(), false));
		//combSet.print();

		//combSet.print("combSet");
		DeskriptorSet<D> result = combSet.getConstrainingProjection(retainedDimensions);

		tcDualisation->pause();
		//result.print("result");
		//printf("Dualize [dim = %d]: start deskriptor count = %d, end deskriptor count = %d\n", getDimension(), getSize(), result.getSize());
		return result;
	}

	template <class D>
	bool operator == (const DeskriptorSet<D>& a, const DeskriptorSet<D>& b) {
		DeskriptorSet<D> aS = a.getSimplifiedDeskriptorSet();
		DeskriptorSet<D> bS = b.getSimplifiedDeskriptorSet();
		//aS.print();
		//bS.print();
		//CInteger(aS.getBidirDeskriptorCount()).print();
		//CInteger(bS.getBidirDeskriptorCount()).print();

		if (aS.getSize() != bS.getSize()) return false;
		if (aS.getBidirDeskriptorCount() != bS.getBidirDeskriptorCount()) return false;

		for (int q = 0; q < aS.getBidirDeskriptorCount(); q++) if (aS[q] != bS[q]) return false;
		for (int q = aS.getBidirDeskriptorCount(); q < aS.getSize(); q++) {	if (!bS.nonNegativelyGenerates(aS[q])) { return false; } }
		for (int q = bS.getBidirDeskriptorCount(); q < bS.getSize(); q++) { if (!aS.nonNegativelyGenerates(bS[q])) { return false; } }

		return true;
	}

	/* For sparse vektors:
	template <class D>
	DeskriptorSet<D> DeskriptorSet<D>::getDual() const {
		typedef typename D::Vector Vector;
		typedef typename Vector::DataPair DataPair;
		std::vector<D*> transposedDeskriptors = std::vector<D*>();
		DeskriptorSet<D> result = DeskriptorSet<D>(getDimension() + getSize());
		for (int q = 0; q < getDimension(); q++) {
			D* desc = new D();
			desc->getVector().appendDataPair(q, (typename D::DeskriptorElementType) 1);
			transposedDeskriptors.push_back(desc);
			result.addDeskriptor(desc);
		}

		printf("Transposing: start size = %d, dim = %d\n", getSize(), getDimension());
		for (int d = 0; d < getSize(); d++) { vector<DataPair>& dataPairs = deskriptors[d]->getVector().getDataPairs();
			for (unsigned int i = 0; i < dataPairs.size(); i++) {
				D* desc = transposedDeskriptors[dataPairs[i].dim];
				desc->getVector().appendDataPair(getDimension() + d, dataPairs[i].value);
				ASSERT(desc->getVector().sorted);
			}
		}

		std::vector<int> retainedDimensions = std::vector<int>();
		for (int q = 0; q < getDimension(); q++) retainedDimensions.push_back(q);

		result = result.getConstrainingProjection(retainedDimensions);
		printf("Dualize [dim = %d]: start deskriptor count = %d, end deskriptor count = %d\n", getDimension(), getSize(), result.getSize());
		return result;
	}*/



	template <class D>
	DeskriptorSet<D>& DeskriptorSet<D>::operator = (const DeskriptorSet<D>& other) {
		if (this != &other) {
			cleanup();
			dimension = other.dimension;
			minimal = other.minimal;
			sorted = other.sorted;
			bidirDeskCount = other.bidirDeskCount;
			for (unsigned int q = 0; q < other.deskriptors.size(); q++) addDeskriptor(*other.deskriptors[q]);
		}

		return *this;
	}

	template <class D>
	void DeskriptorSet<D>::addDeskriptor(D* deskriptor) {
		deskriptors.push_back(deskriptor);
		minimal = false;
		sorted = false;
		bidNormalized = false;
	}

	/*template <class R>
	Deskriptor<R> Deskriptor<R>::getTransposedDeskriptor(const DeskriptorSet<Deskriptor<R> >& deskriptors, int dimension) {
		Vector resultV;

		for (int q = 0; q < deskriptors.getSize(); q++) resultV.appendElement(deskriptors[q][dimension]);

		return Deskriptor<R>(resultV, bidir);
	}*/

	template <class R> R Deskriptor<R>::getGCD() const {
		return vector.getGCD();
	}

	template <class D> void DeskriptorSet<D>::normalizeBidirectionalDeskriptors() {
		if (!bidNormalized) {
			actualNormalizeBidirectionalDeskriptors();
			bidNormalized = true;
		}
	}

	template <class D> void DeskriptorSet<D>::actualNormalizeBidirectionalDeskriptors() {
		sort();
		int bidCount = 0;
		while ((bidCount < getSize()) && (*this)[bidCount].isBidirectional()) bidCount++;

		DeskriptorSet<D>& X = *this;

		pivotDims.clear();
		rank = 0;

		//print();
		for (int j = 0; j < getDimension(); j++) {
			bool pivotFound = false;
			R pivot;
			for (int i = rank; i < bidCount; i++) {
				if (X[i][j] != (R) 0) {
					if (!pivotFound) { // Found it now! :)
						swapDeskriptors(rank, i);
						pivotFound = true;
						pivot = X[rank][j];
						pivotDims.push_back(j);

						if (pivot < (R) 0) { // We must ensure it's positive
							X[rank] = -X[rank];
							pivot = -pivot;
						}
					} else {
						CExtendedExtendedGCDResults<R> xxGCD = extendedExtendedGCD(pivot, X[i][j]);

						D dRank = xxGCD.leftBezoutFactor * X[rank] + xxGCD.rightBezoutFactor * X[i];
						D dI    =  xxGCD.leftLCMCofactor * X[rank] -  xxGCD.rightLCMCofactor * X[i];

						X[rank] = dRank;
						X[i] = dI;

						ASSERT(X[i][j] == (R) 0);

						pivot = xxGCD.gcd;
					}
				}
			}

			if (pivotFound) { // Now reduce the rows above the rank-row
				for (int i = rank - 1; i >= 0; i--) {
					R quotient = X[i][j] / pivot;
					X[i] = X[i] - quotient * X[rank];
				}

				rank++;
			}
		}

		// Remove descriptors that have been nullified
		for (int q = rank; q < bidCount; q++) delete deskriptors[q];
		for (int q = rank; q < getSize() - (bidCount - rank); q++) deskriptors[q] = deskriptors[bidCount + (q - rank)];
		deskriptors.resize(getSize() - (bidCount - rank));

		bidirDeskCount = rank;
	}

	/*template <class D> DeskriptorSet<D> DeskriptorSet<D>::getOrthogonalSet() const {
		typedef typename D::DeskriptorElementType R;

		DeskriptorSet<D> result = *this;
		for (int q = 1; q < getSize(); q++) {
			for (int r = 0; r < q; r++) {
				R rr = result[r] % result[r];
				R rq = -(result[r] % result[q]);
				D rrQ = rr * result[q];
				D rqR = rq * result[r];
				result[q] = (rrQ + rqR).getSimplifiedDeskriptor();
			}
		}

		return result;
	}*/

	template <class R> Deskriptor<R> Deskriptor<R>::getProjection(const DeskriptorSet<Deskriptor<R> >& otherSpaceGenerators, bool preOrthogonalize) const {
		Deskriptor<R> result = *this;

		DeskriptorSet<Deskriptor<R> > oth = (preOrthogonalize) ? otherSpaceGenerators.getOrthogonalSet() : otherSpaceGenerators;

		//print();
		//oth.print();
		for (int q = 0; q < oth.getSize(); q++) {
			//oth[q].print();
			R red = -(oth[q] % result);
			//red.print();
			Deskriptor<R> redV = (red * (oth[q]));
			//redV.print();
			R mul = oth[q] % oth[q];
			//mul.print();
			Deskriptor<R> mulV = mul * result;
			//mulV.print();
			result = mulV + Deskriptor<R>(redV.getVector(), result.isBidirectional());
			result = result.getSimplifiedDeskriptor();
			//result.print();
		}

		//(*this).print();
		//result.print();
		//otherSpaceGenerators.print();
		//oth.print();
		for (int q = 0; q < otherSpaceGenerators.getSize(); q++) {
			//result.print();
			//otherSpaceGenerators[q].print();
			//(result % otherSpaceGenerators[q]).print();
			ASSERT_EQ(result % otherSpaceGenerators[q], CInteger(0));
		}
		//result.print();

		return result;
	}

	template <class D> DeskriptorSet<D> DeskriptorSet<D>::getProjection(const DeskriptorSet<D>& otherSpaceGenerators, bool preOrthogonalize) const {
		DeskriptorSet<D> result = DeskriptorSet<D>(getDimension());

		//otherSpaceGenerators.print();
		DeskriptorSet<D> oth = (preOrthogonalize) ? otherSpaceGenerators.getOrthogonalSet() : otherSpaceGenerators;
		//oth.print();

		for (int q = 0; q < getSize(); q++) {
			D projD = (*this)[q].getProjection(oth, false);
			if (projD != D::getZeroDeskriptor(projD.getSpaceDimension(), projD.isBidirectional())) result.addDeskriptor(projD);
		}

		return result;
	}

	template <class D> DeskriptorSet<D> DeskriptorSet<D>::removeLinearGeneratedVertices(const DeskriptorSet<D>& otherSpaceGenerators, bool preOrthogonalize) const {
		DeskriptorSet<D> result = DeskriptorSet<D>(getDimension());

		//otherSpaceGenerators.print();
		DeskriptorSet<D> oth = (preOrthogonalize) ? otherSpaceGenerators.getOrthogonalSet() : otherSpaceGenerators;
		//oth.print();

		for (int q = 0; q < getSize(); q++) {
			D projD = (*this)[q].getProjection(oth, false);
			if (projD != D::getZeroDeskriptor(projD.getSpaceDimension(), projD.isBidirectional())) result.addDeskriptor((*this)[q]);
		}

		return result;
	}

	template <class R> Deskriptor<R> Deskriptor<R>::getSimplifiedDeskriptor() const {
		return Deskriptor<R>(vector.getGCDNormalizedVektor(), bidir);
	}

	template <class D>
	std::string DeskriptorSet<D>::toString() const {
		std::string result = std::string(getTypeName()) + "(" + CInteger(dimension).toString() + ", ";

		for (unsigned int q = 0; q < deskriptors.size(); q++) result += ((q > 0) ? ", " : "") + deskriptors[q]->toString();

		return result + ")";
	}

	template <class D>
	std::string DeskriptorSet<D>::toStringJavaScript() const {
		std::string desS = "";
		if ((deskriptors.size() != 0) && (dimension != 0)) {
			desS += "new DeskriptorSet(" + CInteger(dimension).toString() + ", [";
			for (unsigned int q = 0; q < deskriptors.size(); q++) { if (q > 0) desS += ", ";
				desS += deskriptors[q]->toStringJavaScript();
			}
			desS += "])";
		}

		return desS;
	}

	template <class D>
	std::string DeskriptorSet<D>::toStringHtml() const {
		std::string desS = "";
		if ((deskriptors.size() != 0) && (dimension != 0)) {
			desS += "<table cellspacing=\"0\" cellpadding=\"0\" border=\"1\">";
			//ITTc(vector<D*>, d, deskriptors) desS += "<tr><td>" + Html::wrapInJS((*d)->toStringJavaScript()) + "</td></tr>";
			ITTc(vector<D*>, d, deskriptors) desS += "<tr><td>" + (*d)->toStringHtml() + "</td></tr>";
			desS += "</table>";
		}

		return Html::wrapper("<table><tr><td>" + getTypeNameHtml() + "</td><td>Dim: "
				+ CInteger(dimension).toString() + "</td></tr></table>",
				   desS, std::string("AABBCC"));
	}

	/* Sparse vektors:
	template <class D>
	DeskriptorSet<D> DeskriptorSet<D>::getTransposedDeskriptorSet() const {
		typedef typename D::Vector Vector;
		typedef typename Vector::DataPair DataPair;
		std::vector<D*> transposedDeskriptors = std::vector<D*>(getDimension(), NULL);

		DeskriptorSet<D> result = DeskriptorSet<D>(getSize());
		printf("Transposing: start size = %d, dim = %d\n", getSize(), getDimension());
		for (int d = 0; d < getSize(); d++) { vector<DataPair>& dataPairs = deskriptors[d]->getVector().getDataPairs();
			for (unsigned int i = 0; i < dataPairs.size(); i++) {
				D* desc = transposedDeskriptors[dataPairs[i].dim];
				if (desc == NULL) {
					desc = new D();
					transposedDeskriptors[dataPairs[i].dim] = desc;
					result.addDeskriptor(desc);
					ASSERT(desc->getVector().sorted);
				}
				desc->getVector().appendDataPair(d, dataPairs[i].value);
				ASSERT(desc->getVector().sorted);
			}
		}

		//for (int q = 0; q < result.getSize(); q++) result[q].getVector().sorted = true;

		printf("Transposing: start size = %d, dim = %d -- end size = %d, dim = %d\n", getSize(), getDimension(), result.getSize(), result.getDimension());
		return result;
	}*/

	template <class R> struct DeskriptorRedundanceStruct {
		Deskriptor<R>* d;
		bool certainlyNonRedundant;

		DeskriptorRedundanceStruct(Deskriptor<R>* iD, bool iCertainlyNonRedundant = false)
		: d(iD), certainlyNonRedundant(iCertainlyNonRedundant) { }
	};

	template <class R> bool pDeskriptorRedundanceStructCompare(const DeskriptorRedundanceStruct<R>& vrs1, const DeskriptorRedundanceStruct<R>& vrs2) {
		return (*(vrs1.d) < *(vrs2.d))
		       ? true
		       : ((*(vrs1.d) == *(vrs2.d)) ? (vrs1.certainlyNonRedundant && !vrs2.certainlyNonRedundant) : false);
	}

	template <class R>
	std::list<DeskriptorRedundanceStruct<R>* > retainNotNonNegativelyGeneratedDeskriptors(std::list<DeskriptorRedundanceStruct<R>* > desks) {
		typedef typename std::list<DeskriptorRedundanceStruct<R>* >::iterator listIt;

		std::list<DeskriptorRedundanceStruct<R>* > nonRed = std::list<DeskriptorRedundanceStruct<R>* >();

		while (desks.size() > 0) {
			DeskriptorSet<Deskriptor<R> > othD = DeskriptorSet<Deskriptor<R> >((*desks.begin())->d->getSpaceDimension());
			//list<Deskriptor<R>* > othD;

			listIt np = desks.begin();
			DeskriptorRedundanceStruct<R>* currD = *np;
			//currD->d->print(); printf(":\n");
			if (!currD->certainlyNonRedundant) {
				for (listIt slit = nonRed.begin(); slit != nonRed.end(); ++slit) othD.addDeskriptor(*(*slit)->d);
				for (listIt slit = ++np; slit != desks.end(); ++slit) othD.addDeskriptor(*(*slit)->d);
				//othD.print();
				if (!othD.nonNegativelyGenerates(*(currD->d))) { nonRed.push_back(currD); /*printf("nonred\n");*/ }
				else { /*printf("red\n");*/ }
			} else { nonRed.push_back(currD); /*printf("cnonred\n");*/ }
			desks.erase(desks.begin());
		}

		return nonRed;
	}

	template <class D> void DeskriptorSet<D>::extractImplicitBidirectionalDeskriptors(bool retainOnlyBidirs) {
		typedef typename R::RationalType Rational;

		//print();
		normalizeBidirectionalDeskriptors();
		//print();

		// Include only unidir desks after projecting out vector space generated by bidirs and add positivity constraint for each of these
		D (*ZD)(int dimension, bool bidir) = &(D::getZeroDeskriptor);
		D (*CD)(int dimension, R constant, bool bidir) = &(D::getConstantDeskriptor);
		D (*UD)(int dimension, int unitDim, bool bidir) = &(D::getUnitDeskriptor);

		for (int s = getDimension() - getBidirSubSet().getSize(); s > 0; s--) {
			DeskriptorSet<D> bidDS = getBidirSubSet();
			//bidDS.print("bidDS");
			DeskriptorSet<D> uniDS = getUnidirSubSet().getProjection(bidDS);
			//uniDS.print("uniDS");
			DeskriptorSet<D> transpUniDS = uniDS.getTransposedDeskriptorSet(true);
			//transpUniDS.print("transpUniDS");
			DeskriptorSet<D> transpDS = transpUniDS << transpUniDS;
			//transpDS.print("transpDS");

			DeskriptorSet<D> equations = DeskriptorSet<D>(transpDS.getDimension() + 1);
			for (int q = 0; q < transpDS.getSize(); q++) equations.addDeskriptor(ZD(1, true) << transpDS[q]);

			DeskriptorSet<D> inequalities = DeskriptorSet<D>(transpDS.getDimension() + 1);
			for (int q = 0; q < transpUniDS.getDimension(); q++) {
				inequalities.addDeskriptor(UD(1 + transpDS.getDimension(), 1 + q, false));
				inequalities.addDeskriptor(UD(1 + transpDS.getDimension(), 1 + transpUniDS.getDimension() + q, false));
			}
			inequalities.addDeskriptor(CD(1, (R) -1, false) << CD(transpUniDS.getDimension(), (R) 1, false) << CD(transpUniDS.getDimension(), (R) 0, false));
			inequalities.addDeskriptor(CD(1, (R) -1, false) << CD(transpUniDS.getDimension(), (R) 0, false) << CD(transpUniDS.getDimension(), (R) 1, false));

			//equations.print("equations");
			//inequalities.print("inequalities");
			//domCons.print();

			LpProblem<R> ipp = LpProblem<R>(transpDS.getDimension());
			ipp.setDomainPolyheder(Flat<Hadron<DeskriptorSet<D> > >(Hadron<DeskriptorSet<D> >(C, inequalities)));
			ipp.setDomainEqualities(Flat<Hadron<DeskriptorSet<D> > >(Hadron<DeskriptorSet<D> >(C, equations)));

			Vektor<Rational>* solutionRat = ipp.solve();
			if (solutionRat == NULL) break;

			Vektor<Rational>& solRat = (*solutionRat);
			//solRat.print();

			R denomLCM = (R) 1;
			for (int q = 0; q < solRat.getLength(); q++) denomLCM = lcm(denomLCM, solRat[q].denominator);
			Vektor<R> sol = Vektor<R>();
			for (int q = 0; q < solutionRat->getLength(); q++) sol.appendElement(solRat[q].numerator * (denomLCM / solRat[q].denominator));
			//sol.print();

			//ASSERT(false);
			/*ParmaMIPProblem plpp = ParmaMIPProblem(domCons.getColumnCount() - 1, (ParmaConstraintSystem) domCons, (ParmaLinearExpression) CVector<R>::getUnitVector(domCons.getColumnCount(), 1), Parma_Polyhedra_Library::MINIMIZATION);
			plpp.solve();
			if (!plpp.is_satisfiable()) break;*/

			//solution->print();
			//CVector<R> sol = CVector<R>::fromParmaGenerator(plpp.optimizing_point());
			//ASSERT(sol[0] == (R) 1);
			//sol = sol.getSubVector(1, sol.getLength() - 1);
			//sol.print();
			Deskriptor<R> kerV = ZD(getDimension(), true);
			for (int q = 0; q < uniDS.getSize(); q++) kerV += sol[q] * uniDS[q];
			//kerV.print();
			kerV.getProjection(bidDS);
			addDeskriptor(kerV.getProjection(bidDS).getSimplifiedDeskriptor());
		}

		DeskriptorSet<D> bidDS = getBidirSubSet();
		//*this = bidDS >> getUnidirSubSet().getProjection(bidDS); // XXX: Not necessary if no implicit bidir descs are found
		if (retainOnlyBidirs) { *this = bidDS; }
		else { *this = bidDS >> getUnidirSubSet().removeLinearGeneratedVertices(bidDS); } // XXX: Not necessary if no implicit bidir descs are found
		//resultCons.print();

		//return CFlat<CModule<R> >(1, CModule<R>(C, matrixToDescriptorSet(resultCons)));
	}

	template <class D>
	DeskriptorSet<D> DeskriptorSet<D>::getSimplifiedDeskriptorSet(int numberOfCertainlyNonRedundant) const {
		std::vector<DeskriptorRedundanceStruct<R> > desks = std::vector<DeskriptorRedundanceStruct<R> >();

		//print();
		DeskriptorSet<D> Y = DeskriptorSet<D>(getDimension());
		for (int q = 0; q < getSize(); q++) {
			Y.addDeskriptor((*this)[q].getSimplifiedDeskriptor());
		}
		Y.sort();
		//Y.print();
		Y.extractImplicitBidirectionalDeskriptors();
		//Y.print();
		Y.normalizeBidirectionalDeskriptors();
		//Y.print();

		// Fill Y with simplified (div by GCD) descriptors
		for (int q = 0; q < Y.getSize(); q++) {
			desks.push_back(DeskriptorRedundanceStruct<R>(&Y[q], false /*q < numberOfCertainlyNonRedundant*/));
		}

		//Y.print();
		// Build a list of less redundant deskriptors (removes duplicates)
		std::list<DeskriptorRedundanceStruct<R>* > remDesks = std::list<DeskriptorRedundanceStruct<R>* >();
		std::vector<DeskriptorRedundanceStruct<R> > desks2 = std::vector<DeskriptorRedundanceStruct<R> >();
		std::sort(desks.begin(), desks.end(), pDeskriptorRedundanceStructCompare<R>);
		D* lastD = NULL;
		DeskriptorSet<D> B = DeskriptorSet<D>(getDimension());
		int unidirs = 0;
		for (unsigned int q = 0; q < desks.size(); q++) {
			if ((lastD == NULL) || (*(desks[q].d) != *lastD)) {
				if (desks[q].d->isBidirectional()) { /*B.addDeskriptor(*desks[q].d);*/ desks2.push_back(desks[q]); }
				else { //desks[q].d->print();
					 desks2.push_back(desks[q]);
					 unidirs++;
				}
				//remDesks.push_back(&desks[q]);
				lastD = desks[q].d;
			}
		}

		//B.print();
		B.bidirDeskCount = 0;

		//for (int q = 0; q < B.getSize(); q++)
		ITT(std::vector<DeskriptorRedundanceStruct<R> >, lit, desks2) {
			remDesks.push_back(&(*lit));
		}
		if (unidirs > 0) remDesks = retainNotNonNegativelyGeneratedDeskriptors(remDesks);

		// XXX: Must re-sort?  Want to change name of simplification to normalisation?
		D* prevD = NULL; // for asserting sortedness
		ITT(std::list<DeskriptorRedundanceStruct<R>* >, lit, remDesks) {
			B.addDeskriptor(*((*lit)->d));
			if ((*lit)->d->isBidirectional()) B.bidirDeskCount = B.bidirDeskCount + 1;
			ASSERT((prevD == NULL) || (*prevD < *(*lit)->d));
			prevD = (*lit)->d;
		}

		ASSERT(B.bidirDeskCount >= 0);

		return B;
	}

	template <class D>
	DeskriptorSet<D> DeskriptorSet<D>::getTransposedDeskriptorSet(bool bidir) const {
		typedef typename D::Vector Vector;

		DeskriptorSet<D> result = DeskriptorSet<D>(getSize());
		for (int d = 0; d < getDimension(); d++) {
		    Vector v = Vector(getSize());
		    for (int e = 0; e < getSize(); e++) v[e] = (*this)[e][d];
		    result.addDeskriptor(D(v, bidir));
		}

		return result;
	}

	template <class D>
	void DeskriptorSet<D>::actualMinimisation() {
		typedef typename DeskriptorType::DeskriptorElementType R;
		if (false) {
			// this one only works for (bidir? --> need glb + ? otherwise) hypercubes
			rank = 0;
			pivotDims = std::vector<int>();
			nonPivotDims = std::vector<int>();
			DeskriptorSet<D>& X = *this;

			for (int j = 0; j < this->getDimension(); j++) {
				if ((j & 0x3F) == 0) printf("%d\n", j);
				bool pivotFound = false;
				for (int i = rank; i < this->getSize(); i++) {
					const R& Xij = X[i](j);
					if (Xij != (R) 0) {
						if (!pivotFound) { // Found it now! :)
							printf("Dim %d, pivot: %d\n", j, i);
							this->swapDeskriptors(rank, i);
							pivotFound = true;
							pivotDims.push_back(j);
						} else X[i] += X[rank];
					}
				}

				if (pivotFound) { // Reduce the rows above the rank-row
					for (int i = rank - 1; i >= 0; i--) X[i] += X[rank]; // XXX: only necessary for normalisation?

					rank++;
				} else nonPivotDims.push_back(j);
			}

			/*
			// Sparse minimisation --> only for sparsevektors
			typedef typename CSparseVector<R>::DataPair DataPair;
			std::vector<std::vector<int > > pivotDeskriptors = std::vector<std::vector<int > >(getDimension());
			for (int d = 0; d < getSize(); d++) {
				vector<DataPair>& dpVect = X[d].getVector().getDataPairs();
				for (unsigned int q = 0; q < dpVect.size(); q++) {
					if (dpVect[q].value != (R) 0) pivotDeskriptors[dpVect[q].dim].push_back(d);
				}
			}

			std::vector<bool> descUsedAsPivot = std::vector<bool>(getSize(), false);
			std::vector<int> pivotDescs = std::vector<int>();
			printf("Minimizing %d dimensional set with %d deskriptors\n", getDimension(), getSize());
			for (int j = 0; j < this->getDimension(); j++) {
				if ((j & 0xFFFF) == 0) printf("%d\n", j);
				bool pivotFound = false;
				std::vector<int>& pivDescs = pivotDeskriptors[j];
				int pivDescI = -1;
				D* pivDesc;
				if (pivDescs.size() > 0) { // pivot desc available?
					for (unsigned int d = 0; d < pivDescs.size(); d++) if (!descUsedAsPivot[pivDescs[d]]) { // available pivot desc found
						descUsedAsPivot[pivDescs[d]] = true;
						pivotFound = true;
						pivDescI = pivDescs[d];
						pivotDescs.push_back(pivDescI);
						pivDesc = &X[pivDescI];
						break;
					}
				}

				if (pivotFound) { // Reduce other pivot deskriptors
					for (unsigned int d = 0; d < pivDescs.size(); d++) if (&X[pivDescs[d]] != pivDesc) { // available pivot desc found
						X[pivDescs[d]] += *pivDesc;
						// XXX: those that are already used as pivots don't need reduction if we don't want to normalize (just minimize)
					}
				}
			}
			// Removing identity deskriptors
			std::vector<D*> retainedDescs = std::vector<D*>();
			for (unsigned int p = 0; p < pivotDescs.size(); p++) {
				retainedDescs.push_back(&X[pivotDescs[p]]);
				deskriptors[pivotDescs[p]] = NULL;
			}
			for (int q = getSize() - 1; q >= 0; q--) if (deskriptors[q] != NULL) delete deskriptors[q];
			deskriptors = retainedDescs;
			*/
			//printf("\n");
		} else { *this = getSimplifiedDeskriptorSet(); }
	}

	template <class D>
	bool DeskriptorSet<D>::bidirectionallyGenerates(const D& deskriptor) const {
		typedef typename D::Vector Vector;

		DeskriptorSet<D> transpDS = getTransposedDeskriptorSet(true);

		//DeskriptorSet (*U)(int dim, bool bid) = &(DeskriptorSet::unitDeskriptorSet);
		//DeskriptorSet (*Z)(int rows, int cols, bool bid) = &(DeskriptorSet::zeroDeskriptorSet);

		DeskriptorSet pointMatrixT = DeskriptorSet(1);
		//deskriptor.print();
		for (int q = 0; q < deskriptor.getVector().getLength(); q++) {
			Vector v = Vector(0);
			v.appendElement(deskriptor.getVector()[q]);
			pointMatrixT.addDeskriptor(D(v, true));
		}
		//pointMatrixT.print();

		//int aDim = getDimension();
		//this->print();
		//transpDS.print();
		//DeskriptorSet diophEq =    (Z(aDim, 1, true) << U(aDim, true)    << -transpDS)
		                        //>> (-pointMatrixT << U(aDim, true)    << Z(aDim, transpDS.getDimension(), true));

		//DeskriptorSet diophEq =    ;
		//diophEq.print();
		DeskriptorSet orthM = (pointMatrixT << -transpDS).getDual();
		if (orthM.getSize() == 0) return false;
		//orthM.print();

		typedef typename D::DeskriptorElementType R;
		R hGcd = orthM[0][0];
		for (int q = 1; q < orthM.getSize(); q++) {
			if (hGcd == 1) return true;
			hGcd = gcd(hGcd, orthM[q][0]);
		}
		return (hGcd.getAbs() == 1);
	}

	/*
	 * Determines redundance of deskriptor
	 * Deskriptor is redundant if it (or a multiple of it) is non-negatively generated by the deskriptors in the set
	 */
	template <class D>
	bool DeskriptorSet<D>::nonNegativelyGenerates(const D& deskriptor) const {
		typedef typename D::DeskriptorElementType R;
		typedef typename R::RationalType Rational;
		typedef typename std::list<CVector<R>* >::iterator listIt;

		ASSERT(bidirDeskCount >= 0);

		DeskriptorSet<D> dumm = *this;
		//dumm.normalizeBidirectionalDeskriptors();
		//print(); deskriptor.print();
		//D (*UD)(int dim, int unitDim, bool bid) = &(D::getUnitDeskriptor);
		D (*ZD)(int dim, bool bidir) = &(D::getZeroDeskriptor);

		if (getSize() == 0) return (deskriptor == ZD(deskriptor.getSpaceDimension(), true)) || (deskriptor == ZD(deskriptor.getSpaceDimension(), false));

		DeskriptorSet<D> ds = *this;
		dumm.addDeskriptor(-deskriptor);
		DeskriptorSet<D> transpDS = dumm.getTransposedDeskriptorSet(true);

		//transpDS.print();

/*		LpProblem<R> ipp = LpProblem<R>(dumm.getSize());
		DeskriptorSet<D> inequalities = DeskriptorSet(transpDS.getDimension() + 1);
		ASSERT(bidirDeskCount >= 0);
		for (int q = 0; q < getSize() + 1; q++) if (!dumm.deskriptors[q]->isBidirectional()) inequalities.addDeskriptor(UD(transpDS.getDimension() + 1, q + 1, false));

		D sum = -D::getUnitDeskriptor(1, 0, false);
		for (int q = 1; q < dumm.getSize(); q++) sum = sum << D::getZeroDeskriptor(1, false);
		sum = sum << D::getUnitDeskriptor(1, 0, false);
		inequalities.addDeskriptor(sum);*/

		DeskriptorSet<D> equations = DeskriptorSet<D>::zeroDeskriptorSet(transpDS.getSize(), 1, true) << transpDS;
/*
		//equations.print(); inequalities.print();

		ipp.setDomainPolyheder(Flat<Hadron<DeskriptorSet<D> > >(Hadron<DeskriptorSet<D> >(C, inequalities)));
		ipp.setDomainEqualities(Flat<Hadron<DeskriptorSet<D> > >(Hadron<DeskriptorSet<D> >(C, equations)));

		Vektor<Rational>* solution = ipp.solve();
		bool nonNegativelyGenerates = (solution != NULL);
		//if (nonNegativelyGenerates) solution->print();
*/
		//equations.print();
		// Now try with glpk
		glp_prob *lp = glp_create_prob();
		glp_set_obj_dir(lp, GLP_MIN);
		glp_add_rows(lp, equations.getSize());
		for (int q = 0; q < equations.getSize(); q++) glp_set_row_bnds(lp, q + 1, GLP_FX, 0.0, 0.0);
		glp_add_cols(lp, equations.getDimension() - 1);
		for (int q = 0; q < getSize() + 1; q++) if ((!dumm.deskriptors[q]->isBidirectional()) || (q == getSize())) {
			//printf("colb : %i = %f\n", q + 1, (q < getSize() ? 0.0 : 1.0));
			glp_set_col_bnds(lp, q + 1, GLP_LO, (q < getSize() ? 0.0 : 1.0), 0.0);
		} else glp_set_col_bnds(lp, q + 1, GLP_FR, 0.0, 0.0);
		for (int q = 0; q < getSize() + 1; q++) glp_set_obj_coef(lp, q + 1, 0.0);
		glp_set_obj_coef(lp, getSize() + 1, 1.0);

		int ia[equations.getSize() * (equations.getDimension() - 1) + 1];
		int ja[equations.getSize() * (equations.getDimension() - 1) + 1];
		double ar[equations.getSize() * (equations.getDimension() - 1) + 1];
		int x = 1;
		for (int k = 0; k < equations.getSize(); k++) {
			for (int l = 1; l < equations.getDimension(); l++) {
				ia[x] = k + 1;
				ja[x] = l;
				ar[x] = equations[k].getVector()[l].toInt();
				//printf("ia[%i] = %i   ja[%i] = %i  ar[%i] = %f\n", x, ia[x], x, ja[x], x, ar[x]);
				x++;
			}
		}
		glp_load_matrix(lp, equations.getSize() * (equations.getDimension() - 1), ia, ja, ar);
		glp_smcp parms;
		glp_init_smcp(&parms);
		parms.msg_lev = GLP_MSG_OFF;
		//parms.
		glp_simplex(lp, &parms);
		int status = glp_get_status(lp);
		double z = glp_get_obj_val(lp);
		bool nNG = ((status == GLP_OPT) && (z >= 1.0))
				   || (status == GLP_UNBND);
		//printf("Solution GLPK = %s (%f, %i)\n", std::string(nNG ? "y" : "n").c_str(), z, status);
		//printf("Solution PIP = %s\n", std::string(nonNegativelyGenerates ? "y" : "n").c_str());
		/*for (int q = 0; q < getSize() + 1; q++) {
			double x = glp_get_col_prim(lp, q + 1);
			printf("col result %i = %f\n", q + 1, x);
		}*/
/*		for (int q = 0; q < getSize() + 1; q++) {
			//printf("colb : %i = %f\n", q + 1, (q < getSize() ? 0.0 : 1.0));
			double x = glp_get_col_prim(lp, q + 1);
			printf("col result %i = %f\n", q + 1, x);
		}*/
		glp_delete_prob(lp);

		//if (nonNegativelyGenerates) (*solution).print();
//		ASSERT(nonNegativelyGenerates == nNG);

		//ParmaMIPProblem plpp = ParmaMIPProblem(domConM.getColumnCount() - 1, (ParmaConstraintSystem) domConM, (ParmaLinearExpression) UV(c.getLength() + 1, 1), Parma_Polyhedra_Library::MINIMIZATION);
		//plpp.solve();

		//print();


		//return nonNegativelyGenerates;
		if (nNG) { nNG = nNG && bidirectionallyGenerates(deskriptor); }

		return nNG;


		//if (solution != NULL) delete solution;
		//return nonNegativelyGenerates;
	}
}



#endif
