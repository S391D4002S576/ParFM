#ifndef RATIONAL_H_
#define RATIONAL_H_

#include <string>
using std::string;

#include "Integer.h"

namespace AlgoTrans {
	template <class R> class CRational {
//		friend CInteger operator + (const CInteger& a, const CInteger& b);
//		friend CInteger operator - (const CInteger& a, const CInteger& b);
//		friend CInteger operator * (const CInteger& a, const CInteger& b);
		template <class G>
		friend bool operator == (const CRational<G>& a, const CRational<G>& b);

		template <class G>
		friend bool operator < (const CRational<G>& a, const CRational<G>& b);
		template <class G>
		friend bool operator > (const CRational<G>& a, const CRational<G>& b);
		
		template <class G>
		friend CRational<G> operator / (const G& a, const G& b);
		
		//friend CInteger operator - (const CInteger& a);
		
		//friend CInteger& operator += (CInteger& a, const CInteger b);
		
		//friend CInteger gcd(const CInteger& a, const CInteger& b);
		//friend CExtendedGCDResults<CInteger> extendedGCD(const CInteger& a, const CInteger& b);
		//friend CExtendedExtendedGCDResults<CInteger> extendedExtendedGCD(const CInteger& a, const CInteger& b);		
	protected:
		
	public:
		typedef R IntegerType;
		
		R numerator, denominator;
		CRational(R iNumerator, R iDenominator) : numerator(iNumerator), denominator(iDenominator) {};
		
		CRational(const CRational& r) { *this = r; }
		CRational(const R& i) { numerator = i; denominator = 1; }
		CRational() { numerator = 0; denominator = 0; };
		
		R getCeil() { return (numerator + denominator - R::getOne()) / denominator; }
		R getFloor() { return numerator / denominator; }
		
		CRational& operator = (R r) {
			numerator = r;
			denominator = 1;

			return *this;
		};
		
		CRational& operator = (const CRational& r) {
			if (this != &r) {
				numerator = r.numerator;
				denominator = r.denominator;
			}
			
			return *this;
		};
		
		string toString() const { return numerator.toString() + "/" + denominator.toString(); }
		
		void print() const;
		
		bool isZero() const { return (denominator != 0) && (numerator == 0); };
		bool isOne() const { return (denominator != 0) && (numerator == denominator); };
		
		static CRational getZero() { return CRational(0); }
		static CRational getOne() { return CRational(1); }
		
		string toStringHtml() const { 
			return "<table><tr><td>" + numerator.toString() + "</td><td>"
				+ denominator.toString() + "</td></tr></table>"; 
		};
	};
 
	template <class R>
	bool operator == (const CRational<R>& a, const CRational<R>& b);
	//bool operator != (const CInteger& a, const CInteger& b);
	
	template <class R>
	bool operator < (const CRational<R>& a, const CRational<R>& b);
	template <class R>
	bool operator > (const CRational<R>& a, const CRational<R>& b);
	
	template <class R>
	CRational<R> operator / (const R& a, const R& b);	
}

#include "Rational.cpp"

#endif
