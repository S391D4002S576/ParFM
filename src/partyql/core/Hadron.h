#ifndef HADRON2_H_
#define HADRON2_H_

#include <vector>
#include <list>
#include <algorithm>
#include <string>
using std::min;

#include "Declarations.h"
#include "Descriptor.h"
#include <typeinfo>

#include "HadronDescription.h"

#include "util/Time.h"
#include "../../../cute/cute.h"

namespace AlgoTrans {

	template <class DeskriptorSet> class Hadron {
	public:
		typedef DeskriptorSet DeskriptorSetType;
		typedef typename DeskriptorSet::DeskriptorType DeskriptorType;
	private:
		void cleanup();
		void cleanDeskriptors(Description description);

		virtual void dualize();

		int dimension;
		DeskriptorSet *constraints, *generators;

		DeskriptorSet*& getDeskriptorsP(Description description) { return (description == G) ? generators : constraints; }
	protected:
		Hadron() : dimension(-1), constraints(NULL), generators(NULL) { }; // For copy-construction in derived classes (???)
		Hadron(int iDimension, DeskriptorSet* iConstraints, DeskriptorSet* iGenerators) : dimension(iDimension), constraints(iConstraints), generators(iGenerators) { };
	public:
		Hadron(int iDimension, Description description) : dimension(iDimension), constraints(NULL), generators(NULL) { getDeskriptorsP(description) = new DeskriptorSet(dimension); };
		Hadron(Description description, const DeskriptorSet& deskriptors) : dimension(deskriptors.getDimension()), constraints(NULL), generators(NULL) { getDeskriptorsP(description) = new DeskriptorSet(deskriptors); };
		Hadron(const Hadron& iOriginal) : dimension(-1), constraints(NULL), generators(NULL) { *this = iOriginal; }
		Hadron& operator = (const Hadron& other);

		virtual ~Hadron() { cleanup(); }

		static Hadron source(int iDimension) { return Hadron(iDimension, G); };
		static Hadron universe(int iDimension) { return Hadron(iDimension, C); };

		int getDimension() const { return dimension; };
		int getDimensionality();

		void setDeskriptors(Description description, const DeskriptorSet& deskriptors);
		const DeskriptorSet& getDeskriptors(Description description);

		bool isDescriptionAvailable(Description description) { return getDeskriptorsP(description) != NULL; };

		Hadron getDual() const;
		Hadron getProjection(std::vector<int> retainedDimensions) const;
		Hadron getDimensionPermutation(std::vector<int> sourceDimensions) const;

		bool isSource();

		//CModule<R> getReflexiveTransitiveClosure(int spaceDim) const;

		// A hadron as deskriptor:
		static Hadron getUnitDeskriptor(int dimension, int unitDimension) {
			DeskriptorSet g = DeskriptorSet(dimension);

			g.addDeskriptor(DeskriptorType::getUnitDeskriptor(dimension, unitDimension));

			return Hadron<DeskriptorSet>(G, g);
		};

		static Hadron getZeroDeskriptor(int dimension) { return source(dimension); };

		bool contains(const DeskriptorType& element) const;

		//static CDeskriptor<R> getTransposedDeskriptor(const CDeskriptorSet<CDeskriptor<R> >& deskriptors, int dimension);

		class MeetOperation { public: static Hadron evaluate(const Hadron& a, const Hadron& b) { return a * b; } };
		class JoinOperation { public: static Hadron evaluate(const Hadron& a, const Hadron& b) { return a + b; } };

		// Text:
		string toString(bool forceBothDescriptions = false);
		string toStringHtml(bool forceBothDescriptions = false);
		string toStringJavaScript(bool forceBothDescriptions = false);

		void print(bool forceBothDescriptions = false) { printf("%s\n", toString(forceBothDescriptions).c_str()); };

		static std::string getTypeName() { return (SHORT_TYPENAMES ? "H<" : "Hadron<" ) + DeskriptorSet::getTypeName() + ">"; }
		static std::string getTypeNameHtml() { return Html::wrapType((SHORT_TYPENAMES ? "H" : "Hadron" ), DeskriptorSet::getTypeNameHtml()); }
	};

	template <class DeskriptorSet>
	Hadron<DeskriptorSet> combineDeskriptorSets(Hadron<DeskriptorSet>& a, Hadron<DeskriptorSet>& b, Description description) {
		ASSERTM(a.getDimension() == b.getDimension(), "Dim(a) = " + CInteger(a.getDimension()).toString() + " != " + CInteger(b.getDimension()).toString() + " = Dim(b)");

		//a.getDeskriptors(description).print();
		//b.getDeskriptors(description).print();
		/*if (a.isDescriptionAvailable(description) && b.isDescriptionAvailable(description))*/ {
		  DeskriptorSet resultDS = a.getDeskriptors(description) + b.getDeskriptors(description);
		  //resultDS.print();
		  resultDS.minimize();
		  //resultDS.print();

		  /*
		   * A <= B ^ C <= B ==> A ^ C <= B
		   * A <= B v C <= B ==>
		   */

		  return Hadron<DeskriptorSet>(description, resultDS);
		} /*else if (a.isDescriptionAvailable(!description) && b.isDescriptionAvailable(!description)) {
			typedef typename DeskriptorSet::DeskriptorType D;

		    if (tcFastMeet == NULL) { tcFastMeet = new TimeCollector("Fast Meet"); timeCollectors.push_back(tcFastMeet); }
		    tcFastMeet->resume();
			DeskriptorSet dsA = a.getDeskriptors(!description);
			DeskriptorSet dsB = b.getDeskriptors(!description);

			int dim = dsB.getDimension();
			DeskriptorSet dsR = DeskriptorSet(2*dsB.getDimension());

			for (int q = 0; q < dsA.getSize(); q++) dsR.addDeskriptor(dsA[q] << D::getZeroDeskriptor(dim, dsA[q].isBidirectional()));
			for (int q = 0; q < dsB.getSize(); q++) dsR.addDeskriptor((-dsB[q]) << dsB[q]);

			vector<int> retainedDimensions;
			for (int q = 0; q < dim; q++) retainedDimensions.push_back(dim + q);

			dsR = dsR.getConstrainingProjection(retainedDimensions);

			Hadron<DeskriptorSet> fast = Hadron<DeskriptorSet>(!description, dsR);
		    tcFastMeet->pause();
			//Hadron<DeskriptorSet> slow = Hadron<DeskriptorSet>(description, a.getDeskriptors(description) + b.getDeskriptors(description));

			//ASSERT_EQ(fast, slow);
			return Hadron<DeskriptorSet>(!description, dsR);
		} *//*else {
			  DeskriptorSet resultDS = a.getDeskriptors(description) + b.getDeskriptors(description);
			  resultDS.minimize();

			  return Hadron<DeskriptorSet>(description, resultDS);
		}*/
	}

	template <class DeskriptorSet> Hadron<DeskriptorSet> operator + (Hadron<DeskriptorSet> a, Hadron<DeskriptorSet> b) {
		if (tcMeet == NULL) { tcMeet = new TimeCollector("Meet"); timeCollectors.push_back(tcMeet); }
		tcMeet->resume();

		// DualAlt: Alternatives here?
		Hadron<DeskriptorSet> result = combineDeskriptorSets(a, b, G);

		tcMeet->pause();

		return result;
	}

	template <class DeskriptorSet> Hadron<DeskriptorSet> operator * (Hadron<DeskriptorSet> a, Hadron<DeskriptorSet> b) {
		// DualAlt: Alternatives here?

		if (tcJoin == NULL) { tcJoin = new TimeCollector("Join"); timeCollectors.push_back(tcJoin); }
		tcJoin->resume();

		// DualAlt: Alternatives here?
		Hadron<DeskriptorSet> result = combineDeskriptorSets(a, b, C);

		tcJoin->pause();

		return result;
	}

	// Product space
	template <class DeskriptorSet> Hadron<DeskriptorSet> operator ^ (Hadron<DeskriptorSet> a, Hadron<DeskriptorSet> b) {
	    if (tcConcatOperation == NULL) { tcConcatOperation = new TimeCollector("Operator ^"); timeCollectors.push_back(tcConcatOperation); }
	    tcConcatOperation->resume();

		// DualAlt: Alternatives here... --> alts accounted for now?

		Description desc = (a.isDescriptionAvailable(C) && b.isDescriptionAvailable(C)) ? C : G;

		int dA = a.getDimension();
		int dB = b.getDimension();
		DeskriptorSet dsA = a.getDeskriptors(desc);
		DeskriptorSet dsB = b.getDeskriptors(desc);

		DeskriptorSet resultDS = DeskriptorSet(dA + dB);
		for (int q = 0; q < dsA.getSize(); q++) resultDS.addDeskriptor(dsA[q] << DeskriptorSet::DeskriptorType::getZeroDeskriptor(dB, dsA[q].isBidirectional()));
		for (int q = 0; q < dsB.getSize(); q++) resultDS.addDeskriptor(DeskriptorSet::DeskriptorType::getZeroDeskriptor(dA, dsB[q].isBidirectional()) << dsB[q]);

	    tcConcatOperation->pause();

		return Hadron<DeskriptorSet>(desc, resultDS);
	}

	template <class DeskriptorSet> Hadron<DeskriptorSet> operator ! (const Hadron<DeskriptorSet>& a) { return a.getDual(); }

	/*template <class DeskriptorSet> Hadron<DeskriptorSet> operator - (const Hadron<DeskriptorSet>& a, const Hadron<DeskriptorSet>& b) {
		return a && b.getComplement();
	}*/

	/*template <class DeskriptorSet> bool operator <= (const Hadron<DeskriptorSet>& a, const Hadron<DeskriptorSet>& b) { // Reflexive SubSetOf
		return (a - b).isSource();
	}*/

	template <class DeskriptorSet> bool operator == (Hadron<DeskriptorSet> a, Hadron<DeskriptorSet> b) {
		// Options:
		// - norm(G(a)) = norm(G(b)) ? // Current implementation
		// - norm(C(a)) = norm(C(b)) ?
		// - norm(G(a) + C(b)) = universe == norm(G(a) * C(b)) = source ?
		// - norm(C(a) + G(b)) = universe == norm(C(a) * G(b)) = source ?
		if (a.isDescriptionAvailable(C) && b.isDescriptionAvailable(C)) {
			return (a.getDeskriptors(C) == b.getDeskriptors(C));
		} else {
			return (a.getDeskriptors(G) == b.getDeskriptors(G));
		}
	}

	template <class DeskriptorSet>
	int Hadron<DeskriptorSet>::getDimensionality() {
		DeskriptorSet ds = getDeskriptors(C);
		ds = ds.getBidirSubSet();
			int b = ds.getSize();

			return getDimension() - b;
	}

	template <class DeskriptorSet>
	void Hadron<DeskriptorSet>::cleanDeskriptors(Description description) {
		DeskriptorSet*& deskriptors = getDeskriptorsP(description);

		if (deskriptors != NULL) delete deskriptors;
		deskriptors = NULL;
	}

	template <class DeskriptorSet>
	void Hadron<DeskriptorSet>::cleanup() {
		cleanDeskriptors(G);
		cleanDeskriptors(C);
	}

	template <class DeskriptorSet>
	bool Hadron<DeskriptorSet>::isSource() {
		if (generators != NULL)  {
			return generators->isEmpty();
		} else {
			ASSERT(false); // must think this over
		}
	}

	template <class DeskriptorSet>
	Hadron<DeskriptorSet>& Hadron<DeskriptorSet>::operator = (const Hadron<DeskriptorSet>& other) {
		if (this != &other) {
			cleanup();
			dimension = other.dimension;
			generators = (other.generators == NULL) ? NULL : new DeskriptorSet(*other.generators);
			constraints = (other.constraints == NULL) ? NULL : new DeskriptorSet(*other.constraints);
		}

		return *this;
	}

	template <class DeskriptorSet>
	void Hadron<DeskriptorSet>::setDeskriptors(Description description, const DeskriptorSet& deskriptors) {
		ASSERT(dimension == deskriptors.getDimension());

		if (getDeskriptorsP(description) != &deskriptors) {
			cleanup();
			getDeskriptorsP(description) = new DeskriptorSet(deskriptors);
		}
	}

	template <class DeskriptorSet>
	const DeskriptorSet& Hadron<DeskriptorSet>::getDeskriptors(Description description) {
		if (getDeskriptorsP(description) == NULL) { // Create it if we don't have it yet
			getDeskriptorsP(description) = new DeskriptorSet(!(*getDeskriptorsP(!description)));
		}

		return *getDeskriptorsP(description);
	}

	template <class DeskriptorSet>
	string Hadron<DeskriptorSet>::toString(bool forceBothDescriptions) {
		std::string result = getTypeName();
		if (forceBothDescriptions){
		    const DeskriptorSet& gen = getDeskriptors(G);
		    const DeskriptorSet& con = getDeskriptors(C);
			return result + "[G: " + gen.toString() + " | C: " + con.toString() + "]";
		} else {
			return result + "[G: " + ((generators != NULL) ? generators->toString() : "NULL")
		                 + " | C: " + ((constraints != NULL) ? constraints->toString() : "NULL") + "]";
		}
	}

	template <class DeskriptorSet>
	string Hadron<DeskriptorSet>::toStringHtml(bool forceBothDescriptions) {
		return Html::wrapInJS(toStringJavaScript(forceBothDescriptions));
		/*std::string result = "<table><tr><td colspan=\"2\">" + getTypeNameHtml() + "</td></tr>";
		if (forceBothDescriptions){
		    const DeskriptorSet& gen = getDeskriptors(G);
		    const DeskriptorSet& con = getDeskriptors(C);
		    return result + "<tr><td>G:</td><td>" + gen.toStringHtml()
			        + "</td>" + "</tr><tr>" + "<td>C:</td><td>"
			       + "</td></tr></table>";
		} else {
			return result + "<tr><td>G:</td><td>"
			       + ((generators != NULL) ? generators->toStringHtml() : "NULL")
			       + "</td>" +  "<td>C:</td><td>"
			       + ((constraints != NULL) ? constraints->toStringHtml() : "NULL")
			       + "</td></tr></table>";
		}*/
	}

	template <class DeskriptorSet>
	string Hadron<DeskriptorSet>::toStringJavaScript(bool forceBothDescriptions) {
		std::string result = "new Hadron(";
		if (forceBothDescriptions){
		    const DeskriptorSet& gen = getDeskriptors(G);
		    const DeskriptorSet& con = getDeskriptors(C);
		    result += gen.toStringJavaScript() + "," + con.toStringJavaScript();
		} else {
			if (generators != NULL) { result += generators->toStringJavaScript(); }
			else { result += "null";}
			result += ",";
			if (constraints != NULL) { result += constraints->toStringJavaScript(); }
			else { result += "null";}
		}
		return result + ")";
	}

	template <class DeskriptorSet> void Hadron<DeskriptorSet>::dualize() {
		DeskriptorSet* qt = generators;
		generators = constraints;
		constraints = qt;
	}

	template <class DeskriptorSet> Hadron<DeskriptorSet> Hadron<DeskriptorSet>::getDual() const {
		Hadron<DeskriptorSet> result = *this;

		result.dualize();

		return result;
	}

	template <class DeskriptorSet>
	Hadron<DeskriptorSet> Hadron<DeskriptorSet>::getDimensionPermutation(std::vector<int> sourceDimensions) const {
		DeskriptorSet* con = (constraints == NULL) ? ((DeskriptorSet*) NULL) : new DeskriptorSet(constraints->getDimensionPermutation(sourceDimensions));
		DeskriptorSet* gen = (generators == NULL)  ? ((DeskriptorSet*) NULL) : new DeskriptorSet(generators->getDimensionPermutation(sourceDimensions));
		return Hadron<DeskriptorSet>(getDimension(), con, gen);
	}

	template <class DeskriptorSet>
	Hadron<DeskriptorSet> Hadron<DeskriptorSet>::getProjection(std::vector<int> retainedDimensions) const {
		if (generators != NULL) {
			return Hadron<DeskriptorSet>(G, generators->getGeneratingProjection(retainedDimensions));
		} else {
			return Hadron<DeskriptorSet>(C, constraints->getConstrainingProjection(retainedDimensions));
		}
	}

	template <class DeskriptorSet>
	bool Hadron<DeskriptorSet>::contains(const DeskriptorType& element) const {
		if (constraints != NULL) {
			for (int q = 0; q < constraints->getSize(); q++) {
				DeskriptorType& c = (*constraints)[q];
				if (c % element < 0) return false;
				if (c.isBidirectional() && (c % element > 0)) return false;
			}

			return true;
		} else { ASSERT(generators != NULL); // generators
			return generators->nonNegativelyGenerates(element);
		}
	}
}

#endif
