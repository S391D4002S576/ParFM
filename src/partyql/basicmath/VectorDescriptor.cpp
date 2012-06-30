#ifdef VECTORDESCRIPTOR_H_

#include "Cone.h"

namespace AlgoTrans {
 	template <class R>
 	CDescriptorSet<CVectorDescriptor<R> > matrixToDescriptorSet(const CMatrix<R>& matrix) {
		CDescriptorSet<CVectorDescriptor<R> > result = CDescriptorSet<CVectorDescriptor<R> >(matrix.getColumnCount());

		for (int q = 0; q < matrix.getRowCount(); q++) {
			result.addDescriptor(CVectorDescriptor<R>(matrix[q]));
		}

		return result;
	}

 	template <class R>
 	CMatrix<R> DescriptorSetToMatrix(const CDescriptorSet<CVectorDescriptor<R> >& descriptors) {
		CMatrix<R> result = CMatrix<R>(descriptors.getDimension());

		for (int q = 0; q < descriptors.getSize(); q++) {
			result.addRow(descriptors[q].vector);
		}

		return result;
	}

	template <class R>
	CUnidirectionalVectorDescriptorSet<R> operator !(CUnidirectionalVectorDescriptorSet<R> a) {
		return (CCone<R>(C, a).getDualizedCone()).getDescriptors(C);
	}

	template <class R>
	void CBidirectionalVectorDescriptorSet<R>::minimize() {
		*this = CBidirectionalVectorDescriptorSet(matrixToDescriptorSet(DescriptorSetToMatrix(*this).getHermiteNormalForm()));
	}

	template <class R>
	void CUnidirectionalVectorDescriptorSet<R>::minimize() {
		ASSERT(false); // separate for cone and polyhedron?  (so what about flats and modules?)
		//*this = matrixToDescriptorSet(DescriptorSetToMatrix(*this).getHermiteNormalForm());
	}
 }

#endif
