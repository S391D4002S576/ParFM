#ifndef HERMITENORMALFORMREDUCTION_H_
#define HERMITENORMALFORMREDUCTION_H_

#include <list>
using std::list;

#include <vector>
using std::vector;

#include "../scalar/ExtendedGCD.h"
#include "HermiteNormalFormReduction.h"

namespace AlgoTrans {
	template <class R> class CMatrix;

	template <class R> class CHermiteNormalFormReducer {
	public:
		int rank;
		CMatrix<R>& X;
		vector<int> pivotColumns;

		CHermiteNormalFormReducer(CMatrix<R>& iMatrix);

		void reduce();

		CMatrix<R>* reduceAndProduceUnitaryMatrix();
	};
}

#include "HermiteNormalFormReduction.cpp"

#endif /* HERMITENORMALFORMREDUCTION_H_ */
