#ifdef RATIONAL_H_ 
#ifndef RATIONAL_CPP_
#define RATIONAL_CPP_

#include "Rational.h"

namespace AlgoTrans {
	template <class R>
	bool operator == (const CRational<R>& a, const CRational<R>& b) {
		return (a.numerator * b.denominator - b.numerator * a.denominator == 0);
	}
//bool operator != (const CInteger& a, const CInteger& b);

	template <class R>
	bool operator < (const CRational<R>& a, const CRational<R>& b) {
		return (a.numerator * b.denominator - b.numerator * a.denominator < 0);
	}
	
	template <class R>
	bool operator > (const CRational<R>& a, const CRational<R>& b) {
		return (a.numerator * b.denominator - b.numerator * a.denominator > 0);
	}

	template <class R>
	CRational<R> operator / (const R& numerator, const R& denominator) {
		return CRational<R>(numerator, denominator);
	}

}

#endif
#endif
