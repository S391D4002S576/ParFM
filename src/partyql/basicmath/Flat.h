

#ifndef FLAT_H_
#define FLAT_H_

#include <string>
using std::string;

#include <vector>
using std::vector;

#include "Hadron.h"
//#include "Descriptors.h"
#include "../../utils/Html.h"

namespace AlgoTrans {
	template <class S> class CFlat: public S {
		template <class T> friend bool operator == (CFlat<T> a, CFlat<T> b);
		template <class T> friend bool isSignatureCompatible(const CFlat<T>& a, const CFlat<T>& b);
	protected:
		int affineness;
		CFlat() : affineness(-1) { }; // For copy-construction in derived classes
	public:
		// Doesn't adjust the SuperSpace accordingly!
		//void setAffineness(int newAffineness) { affineness = newAffineness; S::setSpaceDimension(getSpaceDimension() + affineness); };
		//void setSpaceDimension(int newSpaceDimension) { S::setSpaceDimension(newSpaceDimension + affineness); }

		//CFlat() { };
		//CFlat(int iAffineness) : affineness(iAffineness) { };
		CFlat(int iAffineness, const S& homogeneousSet);
		CFlat(const CFlat<S>& iOriginal) : S(iOriginal) { affineness = iOriginal.getAffineness(); };

		static CFlat<S> source(int iAffineness, int iSpaceDimension);
		static CFlat<S> universe(int iAffineness, int iSpaceDimension);

		virtual ~CFlat() { };

		CFlat<S> getRestriction(vector<int> retainedDimensions) const;
		S getReflexiveTransitiveClosure();

		// Data
		int getAffineness() const { return affineness; };
		int getSpaceDimension() const { return S::getDimension() - getAffineness(); };

		void copySignature(const CFlat<S>& other);

		bool isEmpty();

		// Operators
		CFlat& operator = (const CFlat& other);

		// Debugging
		string toString(bool forceBothDescriptions = false);
		string toStringHtml(bool forceBothDescriptions = false);
		void print(bool forceBothDescriptions = false);
	};


	/*template <> template <class D>
	class CFlat<CHadron<D> > : public CHadron<D>  {
		typedef typename D::DescriptorSet DescriptorSet;
	public:
		CFlat(int iAffineness, HadronDescription description, const DescriptorSet& descriptors)
		: CHadron<D>(description, descriptors), affineness(iAffineness) { }
	};*/

	template <class S>
	bool CFlat<S>::isEmpty() {
		//print();
		S homogeneousProj = S::getProjection(vector<int>(1, 0)); // Project onto homogeneous dimension
		//homogeneousProj.print();

		// Is 1 in the resulting projection?
		return !(homogeneousProj == (homogeneousProj + S(G, S::DescriptorSetType::unitDescriptorSet(1))));
	}

	template <class S> CFlat<S>& CFlat<S>::operator = (const CFlat<S>& other) {
		if (this != &other) {
			affineness = other.affineness;
			S::operator=(other);
		}

		return *this;
	}

	template <class S> bool isSignatureCompatible(const CFlat<S>& a, const CFlat<S>& b) {
		return (a.affineness == b.affineness);
	}

	template <class S>
	void CFlat<S>::copySignature(const CFlat<S>& other) {
		affineness = other.getAffineness();
	}

	template <class S> bool operator == (CFlat<S> a, CFlat<S> b) {
		return isSignatureCompatible(a, b) && ((S) a == (S) b);
	}

	template <class S> string CFlat<S>::toString(bool forceBothDescriptions) {
		return "Flat<" + std::string(S::getTypeName()) +">(" + CInteger(getAffineness()).toString() + ", " + S::toString(forceBothDescriptions) + ")";
	}

	template <class S> string CFlat<S>::toStringHtml(bool forceBothDescriptions) {
		return Html::wrapper("<table><tr><td>Flat&lt;" + S::getTypeNameHtml() + "&gt;</td><td>"
						+ CInteger(getAffineness()).toString() + "</td></tr></table>",
				             S::toStringHtml(forceBothDescriptions), std::string("BBAACC"));
	}

	template <class S> void CFlat<S>::print(bool forceBothDescriptions) {
		printf("%s\n", toString(forceBothDescriptions).c_str());
	}

	template <class S>
	CFlat<S> CFlat<S>::getRestriction(vector<int> retainedDimensions) const {
		for (int q = 0; q < getAffineness(); q++) {
			retainedDimensions.push_back(getSpaceDimension() + q);
		}

		return CFlat(getAffineness(), S::getRestriction(retainedDimensions));
	}

	template <class S>
	S CFlat<S>::getReflexiveTransitiveClosure() {
		ASSERT((getSpaceDimension() % 2) == 0);
		return S::getReflexiveTransitiveClosure(getSpaceDimension() / 2);
	}

	template <class S>
	CFlat<S>::CFlat(int iAffineness, const S& homogeneousSet) : S(homogeneousSet) {
		affineness = iAffineness;
	}

	template <class S>
	CFlat<S> CFlat<S>::source(int iAffineness, int iSpaceDimension) {
		return CFlat<S>(iAffineness, S(iSpaceDimension + iAffineness, G));
	}

	template <class S>
	CFlat<S> CFlat<S>::universe(int iAffineness, int iSpaceDimension) {
		return CFlat<S>(iAffineness, S(iSpaceDimension + iAffineness, C));
	}

	template <class S>
	CFlat<S> operator + (CFlat<S> a, CFlat<S> b) {
		ASSERT(a.getAffineness() == b.getAffineness());

		return CFlat<S>(a.getAffineness(), (S) a + (S) b);
	}

	template <class S>
	CFlat<S> operator * (CFlat<S> a, CFlat<S> b) {
		ASSERT(a.getAffineness() == b.getAffineness());

		return CFlat<S>(a.getAffineness(), (S) a * (S) b);
	}
}

#endif /* FLAT_H_ */

