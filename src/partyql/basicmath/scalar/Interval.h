#ifndef INTERVAL_H_
#define INTERVAL_H_

namespace AlgoTrans {
	template <class R> class CInterval {
		template <class G> friend bool operator == (const CInterval<G> a, const CInterval<G> b);
	private:
	protected:
	public:
		R lowerBound, upperBound;
		bool lowerBounded, upperBounded;
		
		CInterval(bool iLowerBounded, R iLowerBound, bool iUpperBounded, R iUpperBound) {
			lowerBounded = iLowerBounded; 
			upperBounded = iUpperBounded;
			lowerBound = iLowerBound;
			upperBound = iUpperBound;
		}
		
		CInterval(R iLowerBound, R iUpperBound) {
			lowerBounded = true; upperBounded = true;
			lowerBound = iLowerBound;
			upperBound = iUpperBound;
		}
		
		CInterval() {
			lowerBounded = false; upperBounded = false;
		}
		
		CInterval(bool iLowerBounded, R bound) {
			lowerBounded = iLowerBounded; upperBounded = !iLowerBounded;
			
			if (lowerBounded) { lowerBound = bound; }
			else { upperBound = bound; };
		}
		
		CInterval<typename R::IntegerType> asIntegralInterval() {
			return CInterval<typename R::IntegerType>(
					lowerBounded, lowerBound.getCeil(), upperBounded, upperBound.getFloor()
				);
		}
	};
	
	template <class G>
	bool operator == (const CInterval<G> a, const CInterval<G> b);
}

#include "Interval.cpp"

#endif 
