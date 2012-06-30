#ifndef UNITARYTRANSFORMATIONDESCRIPTOR_H_
#define UNITARYTRANSFORMATIONDESCRIPTOR_H_

#include "Descriptor.h"
#include "Hadron.h"

namespace AlgoTrans {
	template <class R> class CMatrix;

	template <class R>
	class CUnitaryTransformationDescriptor : public CDescriptor {
	private:
		CMatrix<R> transformation;
	public:

		CUnitaryTransformationDescriptor(const CMatrix<R>& iTransformation) : transformation(iTransformation) { }

		const CMatrix<R>& getTransformation() { return transformation; }

		typedef CDescriptorTeam<CUnitaryTransformationDescriptor<R> > DescriptorTeam;

		std::string toString() const { return "UnitaryTransformationDescriptor(" + transformation.toString() + ")"; }
	};

	template <class R>
	CDescriptorTeam<CVectorDescriptor<R> > matrixToDescriptorTeam(const CMatrix<R>& matrix);

	template <class R>
	CMatrix<R> descriptorTeamToMatrix(const CDescriptorTeam<CVectorDescriptor<R> >& descriptors);

}

#endif
