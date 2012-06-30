#ifndef INTEGER_CPP_
#define INTEGER_CPP_

#include "Integer.h"

namespace AlgoTrans {
	bool operator == (const CInteger& a, const CInteger& b) {
		return (mpz_cmp(a.data, b.data) == 0);
	}

	bool operator != (const CInteger& a, const CInteger& b) {
		return !(a == b);
	}

	CInteger operator + (const CInteger& a, const CInteger& b) {
		CInteger result = CInteger();

		mpz_add(result.data, a.data, b.data);

		return result;
	}

	CInteger operator % (const CInteger& a, const CInteger& b) {
		CInteger result = CInteger();

		mpz_mod(result.data, a.data, b.data);

		return result;
	}

	CInteger operator - (const CInteger& a, const CInteger& b) {
		CInteger result = CInteger();

		mpz_sub(result.data, a.data, b.data);

		return result;
	}

	CInteger operator * (const CInteger& a, const CInteger& b) {
		CInteger result;

		mpz_mul(result.data, a.data, b.data);

		return result;
	}

	CInteger operator / (const CInteger& a, const CInteger& b) {
		CInteger result;

		mpz_fdiv_q(result.data, a.data, b.data);

		return result;
	}

	bool operator < (const CInteger& a, const CInteger& b) {
		return (mpz_cmp(a.data, b.data) < 0);
	}

	bool operator > (const CInteger& a, const CInteger& b) {
		return (mpz_cmp(a.data, b.data) > 0);
	}

	bool operator >= (const CInteger& a, const CInteger& b) {
		return (mpz_cmp(a.data, b.data) >= 0);
	}

	bool sameSign(const CInteger& a, const CInteger& b) {
		bool aPos = (a > 0);
		bool bPos = (b > 0);

		return (aPos && bPos) || (!aPos && !bPos);
	}

	CInteger operator - (const CInteger& a) {
		CInteger result;

		mpz_neg(result.data, a.data);

		return result;
	}

	CInteger& operator += (CInteger& a, const CInteger& b) {
		mpz_add(a.data, a.data, b.data);

		return a;
	}

	CInteger& operator -= (CInteger& a, const CInteger& b) {
		mpz_sub(a.data, a.data, b.data);

		return a;
	}

	CInteger gcd(const CInteger& a, const CInteger& b) {
		CInteger result;

		mpz_gcd(result.data, a.data, b.data);

		return result;
	}

	CInteger lcm(const CInteger& a, const CInteger& b) {
		CInteger result;

		mpz_lcm(result.data, a.data, b.data);

		return result;
	}

	CInteger CInteger::getLCMCofactor(CInteger& other) const {
		return (other / gcd(*this, other));
	}

	CExtendedGCDResults<CInteger> extendedGCD(const CInteger& a, const CInteger& b) {
		CExtendedGCDResults<CInteger> result;

		mpz_gcdext(result.gcd.data, result.leftBezoutFactor.data, result.rightBezoutFactor.data,
				   a.data, b.data);

		return result;
	}

	CExtendedExtendedGCDResults<CInteger> extendedExtendedGCD(const CInteger& a, const CInteger& b) {
		CExtendedExtendedGCDResults<CInteger> result;

		CExtendedGCDResults<CInteger> xGCD = extendedGCD(a, b);

		mpz_gcdext(result.gcd.data, result.leftBezoutFactor.data, result.rightBezoutFactor.data,
				   a.data, b.data);

		result.leftLCMCofactor = b / result.gcd;
		result.rightLCMCofactor = a / result.gcd;

		return result;
	}

	bool CInteger::isZero() const { return (*this == 0); }

	bool CInteger::isOne() const { return (*this == 1); }

	void CInteger::print() const { printf("%s\n", toString().c_str()); }

}

#endif /* INTEGER_CPP_ */
