#ifdef COMPUTATION_H_
#ifndef VARIABLE_CPP_
#define VARIABLE_CPP_

namespace AlgoTrans {
	template <class R> CVariable<R>::CVariable(CComputation<R>& iComputation)
	: computation(iComputation) {
		computation.registerVariable(*this);
	}

	template <class R> void CVariable<R>::registerReference(CReference<R>& reference) {
		reference.indexInVariable = references.size();

		references.push_back(&reference);
	}

	template <class R> const CVariable<R>& CComputation<R>::getVariable(int index) const {
		return *variables[index];
	}

	template <class R> CVariable<R>& CComputation<R>::addVariable() {
		return *(new CVariable<R>(*this));
	}


}

#endif
#endif
