#ifndef DOMAIN_H_
#define DOMAIN_H_

#include <vector>
using std::vector;

#include <list>
#include <string>

#include "Hadron.h"
#include "Descriptor.h"
#include "Declarations.h"

#include "util/Time.h"
#include "../../../cute/cute.h"

namespace AlgoTrans {
    template <class S> class CDisjunctiveDomain;

	template <class S> class Flat: public S {
	private:
	protected:
		Flat() { }; // For copy-construction in derived classes (XXX: Why not use the other constructors?)
	public:
		typedef S HadronType;

		Flat(const S& homogeneousSet) : S(homogeneousSet) { };
		Flat(const Flat<S>& iOriginal) : S(iOriginal) { };
		Flat(Description description, DeskriptorSet<typename S::DeskriptorType> deskriptors) : S(description, deskriptors) { };

		int getAffineness() const { return 1; }
		int getSpaceDimension() const { return S::getDimension() - getAffineness(); };

		static Flat<S> source(int iSpaceDimension) { return Flat<S>(S(iSpaceDimension + 1, G)); };
		static Flat<S> universe(int iSpaceDimension) { return Flat<S>(S(iSpaceDimension + 1, C)); };

		virtual ~Flat() { };

		Flat<S> getRestriction(vector<int> retainedDimensions) const;
		S getReflexiveTransitiveClosure();

		int getDimensionality();
		Flat<S> getBidirectionalHull();
		Flat<S> getBidirectionalHull_Lim();

		Flat<S> getProjection(std::vector<int> retainedDimensions) const;

		CDisjunctiveDomain<Flat<S> > getComplement();

		bool isEmpty();
		bool contains(const typename S::DeskriptorType& element);

		class MeetOperation: public S::MeetOperation { public: static Flat evaluate(const Flat& a, const Flat& b) { return a * b; } };
		class JoinOperation: public S::JoinOperation { public: static Flat evaluate(const Flat& a, const Flat& b) { return a + b; } };

		// Debugging
		std::string toString(bool forceBothDescriptions = false);
		std::string toStringHtml(bool forceBothDescriptions = false);
		std::string toStringJavaScript(bool forceBothDescriptions = false);
		//std::string toStringJavaScript(bool forceBothDescriptions = false);
		static std::string getTypeName() { return (SHORT_TYPENAMES ? "F<" : "Flat<" ) + S::getTypeName() + ">"; }
		static std::string getTypeNameHtml() { return Html::wrapType((SHORT_TYPENAMES ? "F" : "Flat" ), S::getTypeNameHtml()); }
		void print(bool forceBothDescriptions = false) { printf("%s\n", toString(forceBothDescriptions).c_str()); }
	};

	/*
	 * W must have - getSpaceDimension()
	 *             - getNegative() --> CDisjunctiveDomain<W>
	 *             - Conversion (of singleton) to CDisjunctiveDomain<W>
	 *             - getTypeName()
	 *             - isEmpty()
	 */
	template <class W> class CDisjunctiveDomain {
		//template <class G> friend bool operator == (const CPolyhedralDomain<G>& a, const CPolyhedralDomain<G>& b);
		//template <class G> friend bool operator == (const CPolyheder<G>& a, const CPolyheder<G>& b);
	private:
		int spaceDimension;
		vector<W*> elements;

		void deleteElements() {
			ITT(vector<W*>, e, elements) delete *e;
			elements.clear();
		};

	protected:
		CDisjunctiveDomain() { } // private because we really need the space dimension for now
	public:
		typedef W ElementDataType;

		CDisjunctiveDomain(int iSpaceDimension) : spaceDimension(iSpaceDimension) { };
		CDisjunctiveDomain(const CDisjunctiveDomain& iOriginal);

		CDisjunctiveDomain(const W& iSingleInitialElement) : spaceDimension(iSingleInitialElement.getSpaceDimension()) { addElement(iSingleInitialElement); };

		static CDisjunctiveDomain universe(int iSpaceDimension) { return CDisjunctiveDomain(W::universe(iSpaceDimension)); }
		static CDisjunctiveDomain source(int iSpaceDimension) { return CDisjunctiveDomain(iSpaceDimension); }

		CDisjunctiveDomain getProjection(std::vector<int> retainedDimensions) const;

		void addElement(W* element);
		void removeEmpties();
		void addElement(const W& element) { addElement(new W(element)); }
		const W& getElement(int index) const { return *elements[index]; }
		const W& operator [] (int index) const { return getElement(index); }
		W& operator [] (int index) { return *elements[index]; }

		int getSpaceDimension() const { return spaceDimension; };
		int getElementCount() const { return elements.size(); };

		int getDimensionality() const;

		CDisjunctiveDomain& operator = (const CDisjunctiveDomain& other);

		W getUnidirectionalHull();
		W getBidirectionalHull();

		//operator string();

		bool isEmpty();
		bool isSpaceCompatibleWith(const CDisjunctiveDomain<W>& other) const;

		std::string toStringHtml(bool forceBothDescriptions = false);
		std::string toStringJavaScript(bool forceBothDescriptions = false);
		std::string toString(bool forceBothDescriptions = false);
		static std::string getTypeName() { return (SHORT_TYPENAMES ? "Dom<" : "DisjunctiveDomain<" ) + W::getTypeName() + ">"; }
		static std::string getTypeNameHtml() { return Html::wrapType((SHORT_TYPENAMES ? "Dom" : "DisjunctiveDomain" ), W::getTypeNameHtml()); }
		void print() { printf((toString() + "\n").c_str()); }

		CDisjunctiveDomain getComplement();

		W getHull() const;
	};

	template <class W> int CDisjunctiveDomain<W>::getDimensionality() const {
		int dimensionality = 0;

		for (int q = 0; q < getElementCount(); q++) {
			int d = elements[q]->getDimensionality();
			dimensionality = (dimensionality < d) ? d : dimensionality;
		}

		return dimensionality;
	}

	template <class W> W CDisjunctiveDomain<W>::getBidirectionalHull() {
		W result = W::source(getSpaceDimension());

		for (int q = 0; q < getElementCount(); q++) {
			//result.print();
			//(*this)[q].print();

			if (tcDDBidHull == NULL) { tcDDBidHull = new TimeCollector("DD Bid Hull"); timeCollectors.push_back(tcDDBidHull); }
			tcDDBidHull->resume();

			W BidHull = (*this)[q].getBidirectionalHull();

			tcDDBidHull->pause();

			if (tcDDSum == NULL) { tcDDSum = new TimeCollector("DD Sum"); timeCollectors.push_back(tcDDSum); }
			tcDDSum->resume();

			result = result + BidHull; // XXX: what about empty poly's?

			tcDDSum->pause();
		}

		//result.print();

		return result;
	}

	template <class W> void CDisjunctiveDomain<W>::removeEmpties() {
		int l = 0;
		for (int q = 0; q < getElementCount(); q++) {
			if (elements[q]->isEmpty()) {
				delete elements[q];
			} else elements[l++] = elements[q];
		}

		elements.resize(l);
	}

	template <class W> W CDisjunctiveDomain<W>::getUnidirectionalHull() {
		W result = W::source(getSpaceDimension());

		//removeEmpties();

		//CInteger(getElementCount()).print();
		if (getElementCount() == 1) return (*this)[0];
		if (getElementCount() > 1)  {
	      result = (*this)[0];
	      for (int q = 1; q < getElementCount(); q++) {
	    	  //result.print();
	    	  result = result + (*this)[q]; // XXX: what about empty poly's?
	      }
		}
		//result.print();

		return result;
	}

	template <class S>
	Flat<S> Flat<S>::getProjection(std::vector<int> retainedDimensions) const {
		return Flat<S>(S::getProjection(retainedDimensions));
	}

	// Product space
	template <class S> Flat<S> operator ^ (Flat<S> a, Flat<S> b) { ASSERT(a.getAffineness() == b.getAffineness());
	    if (tcFlatConcatOperation == NULL) { tcFlatConcatOperation = new TimeCollector("Flat Operator ^"); timeCollectors.push_back(tcFlatConcatOperation); }
	    tcFlatConcatOperation->resume();

		typedef typename S::DeskriptorSetType DeskriptorSet;
		typedef typename DeskriptorSet::DeskriptorType Deskriptor;
		typedef typename Deskriptor::Vector Vektor;

		// DualAlt: Alternatives here... --> alts fixed now?

		Description d = C;//(a.isDescriptionAvailable(G) && b.isDescriptionAvailable(G)) ? G : C;

		int dA = a.getSpaceDimension();
		int dB = b.getSpaceDimension();
		//int dC = a.getAffineness();
		DeskriptorSet dsA = a.getDeskriptors(d);
		DeskriptorSet dsB = b.getDeskriptors(d);

		//a.print(); b.print();
		//dsA.print(); dsB.print();

		std::vector<int> retainedDims = std::vector<int>();
		for (int q = 0; q <= 1 + dA + dB; q++) if (q != 1 + dA) retainedDims.push_back(q);

		Vektor (*ZV)(int dim) = &(Vektor::getZeroVektor);

		DeskriptorSet result = DeskriptorSet(1 + dA + dB);
		for (int q = 0; q < dsA.getSize(); q++) {
			Vektor v = dsA[q].getVector().getSubVektor(0, 1) << dsA[q].getVector().getSubVektor(1, dA) << ZV(dB);
			result.addDeskriptor(Deskriptor(v, dsA[q].isBidirectional()));
		}
		for (int q = 0; q < dsB.getSize(); q++) {
			Vektor v = dsB[q].getVector().getSubVektor(0, 1) << ZV(dA) << dsB[q].getVector().getSubVektor(1, dB);
			result.addDeskriptor(Deskriptor(v, dsB[q].isBidirectional()));
		}

		S sE = S(d, result);
		//sE.print(true);

/*		S s = ((S) a) ^ ((S) b);
		//s.print();
		DeskriptorSet homogEq = DeskriptorSet(2 + dA + dB);
		homogEq.addDeskriptor(Deskriptor::getUnitDeskriptor(1 + dA, 0, true) << -Deskriptor::getUnitDeskriptor(1 + dB, 0, true));

	    if (tcFlatConcatOperationMeet == NULL) { tcFlatConcatOperationMeet = new TimeCollector("Flat Operator ^ Meet"); timeCollectors.push_back(tcFlatConcatOperationMeet); }
	    tcFlatConcatOperationMeet->resume();

		s = s * S(C, homogEq);
		//s.print();

		tcFlatConcatOperationMeet->pause();

		//s.print();


	    if (tcFlatConcatOperationProj == NULL) { tcFlatConcatOperationProj = new TimeCollector("Flat Operator ^ Proj"); timeCollectors.push_back(tcFlatConcatOperationProj); }
	    tcFlatConcatOperationProj->resume();

		s = s.getProjection(retainedDims);

		tcFlatConcatOperationProj->pause();*/
		//s.print();
		/*DeskriptorSet resultDS = DeskriptorSet(1 + dA + dB);
		for (int q = 0; q < dsA.getSize(); q++) {
			resultDS.addDeskriptor(dsA[q] << DeskriptorSet::DeskriptorType::getZeroDeskriptor(dB));
		}
		for (int q = 0; q < dsB.getSize(); q++) {
			resultDS.addDeskriptor(dsB[q].getProjection(0, 1) << DeskriptorSet::DeskriptorType::getZeroDeskriptor(dA) << dsB[q].getProjection(1, dB));
		}*/

		tcFlatConcatOperation->pause();

		//s.print(true);
		//ASSERT_EQ(s, sE);

		return Flat<S>(sE);
	}

	template <class S> int Flat<S>::getDimensionality() {
		if (isEmpty()) { return 0; }
		else { return S::getDimensionality() - 1; }
	}

	template <class S> bool Flat<S>::isEmpty() {
	    if (tcEmpty == NULL) { tcEmpty = new TimeCollector("Empty Operation"); timeCollectors.push_back(tcEmpty); }
	    tcEmpty->resume();
		//double prevT = rtclock();
		//print();
		S homogeneousProj = S::getProjection(vector<int>(1, 0)); // Project onto homogeneous dimension
		//homogeneousProj.print();

		// Is 1 in the resulting projection?
		bool nonEmpty = homogeneousProj.contains(S::DeskriptorType::getUnitDeskriptor(1, 0, false));

	    tcEmpty->pause();

		return !nonEmpty;
	}

	template <class S> Flat<S> Flat<S>::getBidirectionalHull_Lim() {
		typedef typename S::DeskriptorSetType DST;
		    DST genDST = S::getDeskriptors(G);

		    DST resultDST = genDST;
		    for (int q = 0; q < genDST.getSize(); q++) resultDST.addDeskriptor(-genDST[q]);

		    return Flat<S>(S(G, resultDST));
	}

	template <class S> Flat<S> Flat<S>::getBidirectionalHull() {
		//if (isEmpty()) return Flat<S>::source(getSpaceDimension());
		typedef typename S::DeskriptorSetType DST;
		typedef typename DST::DeskriptorType D;
		//D (*UD)(int dim, int unitDim, bool bid) = &(D::getUnitDeskriptor);
		//D (*ZD)(int dim, bool bidir) = &(D::getZeroDeskriptor);
		if (S::isDescriptionAvailable(G)) {
		    DST genDST = S::getDeskriptors(G);
		    DST resultDST = genDST;
		    for (int q = 0; q < genDST.getSize(); q++) resultDST.addDeskriptor(-genDST[q]);

		    return Flat<S>(S(G, resultDST));
		} else {
			DST conDST = S::getDeskriptors(C);
			//conDST.print();
			conDST.extractImplicitBidirectionalDeskriptors(true);
			//conDST.print();
			return Flat<S>(S(C, conDST));
		}
	}

	template <class S> bool Flat<S>::contains(const typename S::DeskriptorType& element) {
		typedef typename S::DeskriptorSetType DeskriptorSet;
		typedef typename S::DeskriptorType Deskriptor;

		DeskriptorSet ds = DeskriptorSet(1 + element.getSpaceDimension());
		ds.addDeskriptor(Deskriptor::getUnitDeskriptor(1, 0) << element);
		S singleray = S(G, ds);

		ASSERT(false); // XXX: not yet implemented

		return false;
		//return (((S) (*this)) == ((S) ((*this) + singleray)));
	}

	template <class S> string Flat<S>::toString(bool forceBothDescriptions) {
		return getTypeName() +"(" + CInteger(getAffineness()).toString() + ", " + S::toString(forceBothDescriptions) + ")";
	}

	template <class S> string Flat<S>::toStringHtml(bool forceBothDescriptions) {
		return Html::wrapInJS(toStringJavaScript(forceBothDescriptions));
		/*return Html::wrapper("<table><tr><td>" + getTypeNameHtml() + "</td><td>"
						+ CInteger(getAffineness()).toString() + "</td></tr></table>",
							 Html::wrapInJS(S::toStringJavaScript(forceBothDescriptions)), std::string("BBAACC"));*/
/*		return Html::wrapper("<table><tr><td>" + getTypeNameHtml() + "</td><td>"
						+ CInteger(getAffineness()).toString() + "</td></tr></table>",
							 S::toStringHtml(forceBothDescriptions), std::string("BBAACC"));*/
	}

	template <class S> string Flat<S>::toStringJavaScript(bool forceBothDescriptions) {
		return "new Flat(" + CInteger(getAffineness()).toString() + ","
				+ S::toStringJavaScript(forceBothDescriptions) + ")";
/*		return Html::wrapper("<table><tr><td>" + getTypeNameHtml() + "</td><td>"
						+ CInteger(getAffineness()).toString() + "</td></tr></table>",
							 S::toStringHtml(forceBothDescriptions), std::string("BBAACC"));*/
	}

	template <class S> Flat<S> Flat<S>::getRestriction(vector<int> retainedDimensions) const {
		for (int q = 0; q < getAffineness(); q++) retainedDimensions.push_back(q);

		return Flat(getAffineness(), S::getRestriction(retainedDimensions));
	}

	template <class S> S Flat<S>::getReflexiveTransitiveClosure() { ASSERT((getSpaceDimension() % 2) == 0);
		return S::getReflexiveTransitiveClosure(getSpaceDimension() / 2);
	}

	template <class S> Flat<S> operator + (Flat<S> a, Flat<S> b) { ASSERT(a.getAffineness() == b.getAffineness());
		return Flat<S>((S) a + (S) b);
	}

	// Product space with inhomogeneous constraint handling
	/*template <class S> Flat<S> operator ^ (Flat<S> a, Flat<S> b) {
		ASSERT(a.getAffineness() == b.getAffineness());

	    S im = (S) a ^ (S) b;

	    typedef typename S::DeskriptorSetType DST;
	    typedef typename S::DeskriptorType DT;

	    DST inhomog = DST(im.getDimension());
	    for (int q = 0; q < a.getAffineness(); q++) {
	    	// (bidir) inhomog equation
	    	DT dt = DT(DT::getUnitDeskriptor(im.getDimension(), q).getVector()
	    			   - DT::getUnitDeskriptor(im.getDimension(), a.getAffineness() + a.getSpaceDimension() + q).getVector());
	    	inhomog.addDeskriptor(dt);
	    	inhomog.addDeskriptor(-dt);
	    }

	    im = im * S(C, inhomog);

	    std::vector<int> retainedDims = std::vector<int>();
	    for (int q = 0; q < a.getAffineness() + a.getSpaceDimension(); q++) retainedDims.push_back(q);
	    for (int q = 0; q < b.getSpaceDimension(); q++) retainedDims.push_back(a.getAffineness() + a.getSpaceDimension() + q);
	    im = im.getProjection(retainedDims);

	    im.print();

		return Flat<S>(im);
	}*/

	template <class S> Flat<S> operator * (Flat<S> a, Flat<S> b) { ASSERT(a.getAffineness() == b.getAffineness());
		return Flat<S>((S) a * (S) b);
	}

	template <class W>
	CDisjunctiveDomain<W> CDisjunctiveDomain<W>::getProjection(std::vector<int> retainedDimensions) const {
		CDisjunctiveDomain<W> result = CDisjunctiveDomain<W>(getSpaceDimension());

		for (int q = 0; q < getElementCount(); q++) result.addElement(getElement(q).getProjection(retainedDimensions));

		return result;
	}

	template <class S> CDisjunctiveDomain<Flat<S> > Flat<S>::getComplement() {
		CDisjunctiveDomain<Flat<S> > result = CDisjunctiveDomain<Flat<S> >(getSpaceDimension());

		typedef typename S::DeskriptorType D;
		typedef typename S::DeskriptorSetType DS;
		DS dsCons = (*this).getDeskriptors(C);

		for (int q = 0; q < dsCons.getSize(); q++) {
			if (dsCons[q].isBidirectional()) {
				DS complConDS = DS(getSpaceDimension() + 1);
				complConDS.addDeskriptor(D(-dsCons[q].getVector() - dsCons[q].getVector().getGCD() * D::getUnitDeskriptor(getSpaceDimension() + 1, 0, false).getVector(), false));
				result.addElement(Flat<S>(S(C, complConDS)));

				complConDS = DS(getSpaceDimension() + 1);
				complConDS.addDeskriptor(D(dsCons[q].getVector() - dsCons[q].getVector().getGCD() * D::getUnitDeskriptor(getSpaceDimension() + 1, 0, false).getVector(), false));
				result.addElement(Flat<S>(S(C, complConDS)));
			} else {
				DS complConDS = DS(getSpaceDimension() + 1);
				complConDS.addDeskriptor(D(-dsCons[q].getVector() - dsCons[q].getVector().getGCD() * D::getUnitDeskriptor(getSpaceDimension() + 1, 0, false).getVector(), false));
				result.addElement(Flat<S>(S(C, complConDS)));
			}
		}
		//result.print();

		return result;
	}

	//template <class W> bool operator == (const CDisjunctiveDomain<W>& a, const CDisjunctiveDomain<W>& b);
	//template <class W> bool operator != (const CDisjunctiveDomain<W>& a, const CDisjunctiveDomain<W>& b);


	template <class W> CDisjunctiveDomain<W>& CDisjunctiveDomain<W>::operator = (const CDisjunctiveDomain<W>& other) {
		if (this != &other) {
			spaceDimension = other.spaceDimension;

			deleteElements();
			for (unsigned int q = 0; q < other.elements.size(); q++) addElement(*other.elements[q]);
		}

		return *this;
	}

	template <class W> bool CDisjunctiveDomain<W>::isEmpty() {
		ITT(vector<W*>, e, elements) if (!((*e)->isEmpty())) return false; // XXX: under some minimisation assumptions, this can be done faster

		return true;
	}

	template <class S> Flat<S> operator && (const Flat<S>& a, const Flat<S>& b) {
		return Flat<S>(((S) a) * ((S) b));
	}

	template <class S> Flat<S> operator || (const Flat<S>& a, const Flat<S>& b) {
		return Flat<S>(((S) a) + ((S) b));
	}

	template <class G> CDisjunctiveDomain<G> operator && (const CDisjunctiveDomain<G>& a, const CDisjunctiveDomain<G>& b) {
		int newSpaceDim = -1;
		if (a.getElementCount() > 0) newSpaceDim = a.getSpaceDimension();
		if (b.getElementCount() > 0) newSpaceDim = b.getSpaceDimension();

		CDisjunctiveDomain<G> result = CDisjunctiveDomain<G>(newSpaceDim);
		for (int q = 0; q < a.getElementCount(); q++) {
			for (int r = 0; r < b.getElementCount(); r++) {
				result.addElement(a[q] && b[r]);
			}
		}

		return result;
	}

	template <class G> CDisjunctiveDomain<G> operator ^ (const CDisjunctiveDomain<G>& a, const CDisjunctiveDomain<G>& b) {
		CDisjunctiveDomain<G> result = CDisjunctiveDomain<G>(a.getSpaceDimension() + b.getSpaceDimension());

		for (int q = 0; q < a.getElementCount(); q++) {
			for (int r = 0; r < b.getElementCount(); r++) result.addElement(a[q] ^ b[r]);
		}

		return result;
	}

	template <class G> CDisjunctiveDomain<G> operator && (const G& a, const CDisjunctiveDomain<G>& b) { return ((CDisjunctiveDomain<G>) a) && b; }
	template <class G> CDisjunctiveDomain<G> operator && (const CDisjunctiveDomain<G>& a, const G& b) { return a && ((CDisjunctiveDomain<G>) b); }

	template <class G> CDisjunctiveDomain<G> operator || (const CDisjunctiveDomain<G>& a, const CDisjunctiveDomain<G>& b) {
		CDisjunctiveDomain<G> result = CDisjunctiveDomain<G>(a.getSpaceDimension());

		for (int q = 0; q < a.getElementCount(); q++) result.addElement(a[q]);
		for (int q = 0; q < b.getElementCount(); q++) result.addElement(b[q]);

		return result;
	}

	template <class W> void CDisjunctiveDomain<W>::addElement(W* element) {
		if (getElementCount() == 0) spaceDimension = element->getSpaceDimension();
		ASSERT(element->getSpaceDimension() == spaceDimension);

		elements.push_back(element);
	}

	template <class W> CDisjunctiveDomain<W> CDisjunctiveDomain<W>::getComplement() {
		CDisjunctiveDomain<W> result = CDisjunctiveDomain<W>(getSpaceDimension());

		result.addElement(W::universe(getSpaceDimension()));

		ITT(vector<W*>, e, elements) result = result && (**e).getComplement();

		return result;
	}

	template <class H> CDisjunctiveDomain<Flat<H> > operator - (Flat<H> a, Flat<H> b) {
		return CDisjunctiveDomain<Flat<H> >(a) && b.getComplement(); // Probably not the most efficient way --> avoid dualisation
	}

	template <class G> CDisjunctiveDomain<G> operator - (CDisjunctiveDomain<G> a, CDisjunctiveDomain<G> b) {
		return a && b.getComplement(); // Probably not the most efficient way --> avoid dualisation
	}

	template <class H> bool operator <= (Flat<H> a, Flat<H> b) { // Reflexive SubSetOf
		//a.print();
		//b.print();
		//(a - b).print();
		return (a - b).isEmpty(); // Probably not the most efficient way
	}

	template <class G> bool operator <= (CDisjunctiveDomain<G> a, CDisjunctiveDomain<G> b) { // Reflexive SubSetOf
		//a.print();
		//b.print();
		//(a - b).print();
		return (a - b).isEmpty(); // Probably not the most efficient way
	}

	template <class H> bool operator == (Flat<H> a, Flat<H> b) {
		//a.print();
		//b.print();
		//if (a <= b) printf("a <= b\n");
		//if (b <= a) printf("b <= a\n");
		if (a.isDescriptionAvailable(C) && b.isDescriptionAvailable(C)) {
			return a.getDeskriptors(C) == b.getDeskriptors(C);
		} else /*if (a.isDescriptionAvailable(G) && b.isDescriptionAvailable(G))*/ {
			//a.getDeskriptors(G).print();
			//b.getDeskriptors(G).print();
			return a.getDeskriptors(G) == b.getDeskriptors(G);
		}

		//return ((H) a == (H) b);
		/*a.getDeskriptors(C).getSimplifiedDeskriptorSet().print();
		b.getDeskriptors(C).getSimplifiedDeskriptorSet().print();
		a.print();
		b.print();
		return a.getDeskriptors(C) == b.getDeskriptors(C); // Faster than previous version?*/
	}

	template <class H> bool operator != (Flat<H> a, Flat<H> b) { return !(a == b); }

	template <class G> bool operator == (CDisjunctiveDomain<G> a, CDisjunctiveDomain<G> b) {
		//a.print();
		//b.print();
		//if (a <= b) printf("a <= b\n");
		//if (a <= b) printf("b <= a\n");
		return (a <= b) && (b <= a); // Probably not the most efficient way --> better: minimize both &  sort&compare
	}

	template <class G> bool operator != (CDisjunctiveDomain<G> a, CDisjunctiveDomain<G> b) { return !(a == b); }

	template <class W> std::string CDisjunctiveDomain<W>::toString(bool forceBothDescriptions) {
		std::string result = getTypeName() + "(";

		for (int q = 0; q < getElementCount(); q++) {
			result += ((q != 0) ? ", " : "") + (*this)[q].toString(forceBothDescriptions);
			ASSERT_EQ(CInteger((*this)[q].getSpaceDimension()), CInteger(getSpaceDimension()));
		}

		return result + ")";
	}

	template <class W> std::string CDisjunctiveDomain<W>::toStringHtml(bool forceBothDescriptions) {
		return Html::wrapInJS(toStringJavaScript(forceBothDescriptions));
		/*std::string result = "<table><tr><td>" + getTypeNameHtml() + "</td><td>";

		for (int q = 0; q < getElementCount(); q++) {
			result += ((q != 0) ? "</td><td>, </td><td>" : "") + Html::wrapInJS((*this)[q].toStringJavaScript(forceBothDescriptions)) + "</td><td>";
		}

		return result + "</td></tr></table>";*/
	}

	template <class W> std::string CDisjunctiveDomain<W>::toStringJavaScript(bool forceBothDescriptions) {
		std::string result = "new DisjunctiveDomain([";

		for (int q = 0; q < getElementCount(); q++) {
			result += (*this)[q].toStringJavaScript(forceBothDescriptions) + ",";
		}

		return result + "])";
	}

	template <class R> bool CDisjunctiveDomain<R>::isSpaceCompatibleWith(const CDisjunctiveDomain<R>& other) const {
		return (getElementCount() == 0) || (other.getElementCount() == 0) || (getSpaceDimension() == other.getSpaceDimension());
	}

	template <class W> W CDisjunctiveDomain<W>::getHull() const {
		W result = W::source(spaceDimension());
		ASSERT(false);

		bool empty = true;
		ITT(vector<W*>, e, elements) if (!((*e)->isEmpty())) {
			if (empty) {
				result = **e;
				empty = false;
			} else {
				result = result + **e;
			}
		} // XXX: else { Remove the element since we know we don't need it? (const...) }

		return result;
	}

	template <class W>
	CDisjunctiveDomain<W>::CDisjunctiveDomain(const CDisjunctiveDomain& iOriginal) {
		*this = iOriginal;

		for (unsigned int q = 0; q < elements.size(); q++) {
			elements[q] = new W(*elements[q]);
		}
	};
}

#endif /* POLYHEDER_H_ */
