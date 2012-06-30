#ifdef INTERVAL_H_
#ifndef INTERVAL_CPP_
#define INTERVAL_CPP_

#include "Interval.h"

namespace AlgoTrans {
	template <class G>
	bool operator == (const CInterval<G> a, const CInterval<G> b) {
		return (a.lowerBounded == b.lowerBounded)
			&& (a.upperBounded == b.upperBounded)
			&& (a.lowerBound == b.lowerBound)
			&& (a.upperBound == b.upperBound);
	}
	
	template <class G>
	CInterval<G> cover(const CInterval<G> a, const CInterval<G> b) {
		CInterval<G> result;
		
		result.lowerBounded = a.lowerBounded && b.lowerBounded;
		result.upperBounded = a.upperBounded && b.upperBounded;
		
		if (result.lowerBounded) result.lowerBound = min(a.lowerBound, b.lowerBound);
		if (result.upperBounded) result.upperBound = max(a.upperBound, b.upperBound);
		
		return result;
	}
}

#endif
#endif
