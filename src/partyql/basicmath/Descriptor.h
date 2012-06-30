#ifndef DESCRIPTOR_H_
#define DESCRIPTOR_H_

#include <vector>
using std::vector;

#include <list>
#include <string>
using std::list;

#include <algorithm>
using std::min;

#include "scalar/Integer.h"
#include "../../utils/Html.h"

namespace AlgoTrans {
	class CDescriptor {
		//typedef DescriptorSet* DescriptorSetP;
	public:
		virtual std::string toString() const = 0;
		void print() const { printf("%s\n", toString().c_str()); };

		// CDescriptor getUnitDescriptor(int dimension, int unitDimension) = 0;

		virtual ~CDescriptor() { };
	};

	template <class D>
	class CDescriptorSet {
	private:
		int dimension;
		vector<D* > descriptors;
		void cleanup();

		bool minimal;
		virtual void actualMinimisation() { ASSERT(false); /* convenient abstract virtuality */ };
	public:
		CDescriptorSet(int iDimension);
		CDescriptorSet(const CDescriptorSet& iOriginal);
		CDescriptorSet& operator = (const CDescriptorSet& other);

		~CDescriptorSet() { cleanup(); };

		const D& operator [] (int i) { return *descriptors[i]; };
		const D& operator [] (int i) const { return *descriptors[i]; };

		int getSize() const { return descriptors.size(); }
		int getDimension() const { return dimension; };

		void addDescriptor(const D& descriptor); // Ownership of descriptor remains with caller
		void addDescriptor(D* descriptor); // Ownership of descriptor transferred to team

		//CDescriptorSet getDual() const;

		void minimize();

		static CDescriptorSet<D> unitDescriptorSet(int dimension);

		std::string toString() const;
		std::string toStringHtml() const;
		void print() const { printf("%s\n", toString().c_str()); };

		static std::string getTypeName() { return "CDescriptorSet<" + D::getTypeName() + ">"; }
		static std::string getTypeNameHtml() { return "CDescriptorSet&lt;" + D::getTypeNameHtml() + "&gt;"; }
	};

	template <class D>
	CDescriptorSet<D> operator + (const CDescriptorSet<D>& a, const CDescriptorSet<D>& b) {
		ASSERT(a.getDimension() == b.getDimension());

		CDescriptorSet<D> result = CDescriptorSet<D>(a.getDimension());

		for (int q = 0; q < a.getSize(); q++) result.addDescriptor(a[q]);
		for (int q = 0; q < b.getSize(); q++) result.addDescriptor(b[q]);

//		result.minimize();

		return result;
	}

	template <class D>
	CDescriptorSet<D> CDescriptorSet<D>::unitDescriptorSet(int dimension) {
		CDescriptorSet<D> result = CDescriptorSet<D>(dimension);

		for (int q = 0; q < dimension; q++) result.addDescriptor(D::getUnitDescriptor(dimension, q));

		return result;
	}

	template <class D>
	void CDescriptorSet<D>::cleanup() {
		for (unsigned int q = 0; q < descriptors.size(); q++) delete descriptors[q];
		descriptors.clear();
	}

	template <class D>
	CDescriptorSet<D>::CDescriptorSet(int iDimension) : dimension(iDimension), minimal(true) {
	}

	template <class D>
	CDescriptorSet<D>::CDescriptorSet(const CDescriptorSet<D>& iOriginal) {
		*this = iOriginal;
	}

	template <class D>
	void CDescriptorSet<D>::minimize() {
		if (!minimal) {
			actualMinimisation();
			minimal = true;
		}
	}

	template <class D>
	CDescriptorSet<D>& CDescriptorSet<D>::operator = (const CDescriptorSet<D>& other) {
		if (this != &other) {
			cleanup();
			dimension = other.dimension;
			minimal = other.minimal;
			for (unsigned int q = 0; q < other.descriptors.size(); q++) addDescriptor(*other.descriptors[q]);
		}

		return *this;
	}

	template <class D>
	void CDescriptorSet<D>::addDescriptor(const D& descriptor) {
		addDescriptor(new D(descriptor));
	}

	template <class D>
	void CDescriptorSet<D>::addDescriptor(D* descriptor) {
		descriptors.push_back(descriptor);
		minimal = false;
	}

	template <class D>
	std::string CDescriptorSet<D>::toString() const {
		std::string result = "DescriptorSet&lt;" + std::string(D::getTypeName()) + "&gt;(" + CInteger(dimension).toString() + ", ";

		for (unsigned int q = 0; q < descriptors.size(); q++) {
			if (q > 0) result += ", ";
			result += descriptors[q]->toString();
		}

		return result + ")";
	}

	template <class D>
	std::string CDescriptorSet<D>::toStringHtml() const {
		std::string desS = "";
		if ((descriptors.size() != 0) && (dimension != 0)) {
			desS += "<table cellspacing=\"0\" border=\"1\">";
			for (unsigned int r = 0; r < descriptors.size(); r++) {
				desS += "<tr><td>" + descriptors[r]->toStringHtml() + "</td></tr>";
			}
			desS += "</table>";
		}

		return Html::wrapper("<table><tr><td>DescriptorSet&lt;" + D::getTypeNameHtml() + "&gt;</td><td>Dim: "
				+ CInteger(dimension).toString() + "</td></tr></table>",
				   desS, std::string("AABBCC"));
	}
}

#endif
