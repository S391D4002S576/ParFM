#ifndef VECTORDESCRIPTOR_H_
#define VECTORDESCRIPTOR_H_

#include "../cute/cute.h"

#include "Descriptor.h"
#include "Hadron.h"

namespace AlgoTrans {
	template <class R> class CMatrix;

	template <class R>
	class CVectorDescriptor : public CDescriptor {
	public:
		CVector<R> vector;

		CVectorDescriptor(const CVector<R>& iVector) : vector(iVector) { }

		static CVectorDescriptor<R> getUnitDescriptor(int dimension, int unitDimension) { return CVectorDescriptor<R>(CVector<R>::getUnitVector(dimension, unitDimension)); };

		std::string toString() const { return vector.toString(); }
		std::string toStringHtml() const { return vector.toStringHtml(); }

		static std::string getTypeName() { return "CVectorDescriptor<" + R::getTypeName() + ">"; }
		static std::string getTypeNameHtml() { return "CVectorDescriptor&lt;" + R::getTypeNameHtml() + "&gt;"; }
	};

	template <class R>
	class CBidirectionalVectorDescriptorSet : public CDescriptorSet<CVectorDescriptor<R> > {
	public:
		typedef CBidirectionalVectorDescriptorSet<R>* Pointer;

		CBidirectionalVectorDescriptorSet(const CDescriptorSet<CVectorDescriptor<R> >& iOriginal) : CDescriptorSet<CVectorDescriptor<R> >(iOriginal) { };

		virtual void minimize();
	};

	template <class R>
	class CUnidirectionalVectorDescriptorSet : public CDescriptorSet<CVectorDescriptor<R> >{
	public:
		typedef CUnidirectionalVectorDescriptorSet<R>* Pointer;

		CUnidirectionalVectorDescriptorSet(const CDescriptorSet<CVectorDescriptor<R> >& iOriginal) : CDescriptorSet<CVectorDescriptor<R> >(iOriginal) { };

		virtual void minimize();
	};

	template <class R>
	CDescriptorSet<CVectorDescriptor<R> > matrixToDescriptorSet(const CMatrix<R>& matrix);

	template <class R>
	CMatrix<R> DescriptorSetToMatrix(const CDescriptorSet<CVectorDescriptor<R> >& descriptors);

	template <class R>
	CDescriptorSet<CVectorDescriptor<R> > operator !(CDescriptorSet<CVectorDescriptor<R> > a) {
		ASSERT(false); // Depends on directionality, so deferred to derived classes

		return a; /// XXX
	}

	template <class R>
	CBidirectionalVectorDescriptorSet<R> operator !(CBidirectionalVectorDescriptorSet<R> a) {
		return matrixToDescriptorSet(!DescriptorSetToMatrix(a));
	}

	template <class R> CUnidirectionalVectorDescriptorSet<R> operator !(CUnidirectionalVectorDescriptorSet<R> a);

	template <class R>
	CUnidirectionalVectorDescriptorSet<R> operator + (CUnidirectionalVectorDescriptorSet<R>& a, CUnidirectionalVectorDescriptorSet<R>& b) {
		return CUnidirectionalVectorDescriptorSet<R>((CDescriptorSet<CVectorDescriptor<R> >&) a + (CDescriptorSet<CVectorDescriptor<R> >&) b);
	}

	template <class R>
	CBidirectionalVectorDescriptorSet<R> operator + (CBidirectionalVectorDescriptorSet<R>& a, CBidirectionalVectorDescriptorSet<R>& b) {
		return CBidirectionalVectorDescriptorSet<R>((CDescriptorSet<CVectorDescriptor<R> >&) a + (CDescriptorSet<CVectorDescriptor<R>& >) b);
	}
}

#include "VectorDescriptor.cpp"

#endif /*VECTORDESCRIPTOR_H_*/
