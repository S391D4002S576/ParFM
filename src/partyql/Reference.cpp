#ifdef COMPUTATION_H_
#ifndef REFERENCE_CPP
#define REFERENCE_CPP

#include "../cute/cute.h"

#include "RADGraph.h"

namespace AlgoTrans {
	template <class R> CReference<R>::CReference(CComputation<R>& iComputation,
			CStatement<R>& iStatement, CVariable<R>& iVariable, UReadWrite iReadWrite, const CAffineTransformation<R>& iTransformation) 
	: computation(iComputation), statement(iStatement), variable(iVariable), readWrite(iReadWrite), transformation(iTransformation) {
		computation.registerReference(*this);
	}
}

#endif
#endif
