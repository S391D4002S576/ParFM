#ifndef HADRON_H_
#define HADRON_H_

#include <vector>
using std::vector;

#include <list>
using std::list;

#include <algorithm>
#include <string>
using std::min;

#include "HadronDescription.h"
#include <typeinfo>

namespace AlgoTrans {

	template <class DescriptorSet>
	class CHadron {
	public:
		typedef DescriptorSet DescriptorSetType;
		typedef typename DescriptorSet::Pointer DescriptorSetP;
	private:
		typedef HadronDescription Description;
	private:
		void cleanup();
		void cleanDescriptors(Description description);

		virtual void dualize();

		int dimension;
		DescriptorSet* constraints;
		DescriptorSet* generators;

		DescriptorSetP& getDescriptorsP(Description description);
	protected:
		CHadron() : dimension(-1), constraints(NULL), generators(NULL) { }; // For copy-construction in derived classes
	public:
		CHadron(int iDimension, Description description);
		CHadron(Description description, const DescriptorSet& descriptors);

		CHadron(const CHadron& iOriginal);
		CHadron& operator = (const CHadron& other);

		~CHadron();

		//void setDimension(int newDimension);
		int getDimension() const { return dimension; };

		void setDescriptors(Description description, const DescriptorSet& descriptors);
		const DescriptorSet& getDescriptors(Description description);

		static CHadron source(int iDimension) { return CHadron(iDimension, G); };
		static CHadron universe(int iDimension) { return CHadron(iDimension, C); };

		CHadron getDual() const;

		// Projection represented in lower-dimensional (projected) space
		//CModule<R> getRestriction(vector<int> retainedDimensiones) const;
		//CModule<R> getReflexiveTransitiveClosure(int spaceDim) const;

		//int getGeneratingVectorCount() const { return matrix.getRowCount(); };

		string toString(bool forceBothDescriptions = false);
		string toStringHtml(bool forceBothDescriptions = false);

		void print(bool forceBothDescriptions = false) { printf("%s\n", toString(forceBothDescriptions).c_str()); };

		static std::string getTypeName() { return "Hadron<" + DescriptorSet::getTypeName() + ">"; }
		static std::string getTypeNameHtml() { return "Hadron&lt;" + DescriptorSet::getTypeNameHtml() + "&gt;"; }
	};

	template <class DescriptorSet>
	CHadron<DescriptorSet> combineHadrons(CHadron<DescriptorSet>& a, CHadron<DescriptorSet>& b, HadronDescription description) {
		ASSERTM(a.getDimension() == b.getDimension(), "Dim(a) = " + CInteger(a.getDimension()).toString() + " != " + CInteger(b.getDimension()).toString() + " = Dim(b)");

		DescriptorSet resultDS = a.getDescriptors(description) + b.getDescriptors(description);
		//resultDS.print();
		resultDS.minimize();
		//resultDS.print();

		return CHadron<DescriptorSet>(description, resultDS);
	}

	template <class DescriptorSet> CHadron<DescriptorSet> join(CHadron<DescriptorSet>& a, CHadron<DescriptorSet>& b) {
		return combineHadrons(a, b, G);
	}

	template <class DescriptorSet> CHadron<DescriptorSet> operator +(CHadron<DescriptorSet> a, CHadron<DescriptorSet> b) {
		return join(a, b);
	}

	template <class DescriptorSet> CHadron<DescriptorSet> meet(CHadron<DescriptorSet>& a, CHadron<DescriptorSet>& b) {
		return combineHadrons(a, b, C);
	}

	template <class DescriptorSet> CHadron<DescriptorSet> operator *(CHadron<DescriptorSet> a, CHadron<DescriptorSet> b) {
		return meet(a, b);
	}

	template <class DescriptorSet> CHadron<DescriptorSet> operator ! (CHadron<DescriptorSet> a) {
		return a.getDual();
	}

	/* For bidirectional descriptors:
	template <class DescriptorSet> bool operator == (CHadron<DescriptorSet> a, CHadron<DescriptorSet> b) {
		// Options:
		// - norm(G(a)) = norm(G(b)) ? // Current implementation
		// - norm(C(a)) = norm(C(b)) ?
		// - norm(G(a) + C(b)) = universe = norm(G(a) * C(b))
		// - norm(C(a) + G(b)) = universe = norm(C(a) * G(b))
		DescriptorSet aDS = a.getDescriptors(G);
		aDS.minimize();
		DescriptorSet bDS = b.getDescriptors(G);
		bDS.minimize();

		if (aDS.getSize() != bDS.getSize()) return false;

		for (int q = 0; q < aDS.Size(); q++) if (aDS[q] != bDS[q]) return false;

		return true;
	}
	*/
}

#include "Hadron.cpp"

#endif
