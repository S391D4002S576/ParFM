#ifndef INTEGER_H_
#define INTEGER_H_

#include <gmp.h>
#include <gmpxx.h>
#include <string>
#include "ExtendedGCD.h"
using std::string;

namespace AlgoTrans {
	template <class R> class CRational;

	class CInteger  {
		friend CInteger operator + (const CInteger& a, const CInteger& b);
		friend CInteger operator - (const CInteger& a, const CInteger& b);
		friend CInteger operator * (const CInteger& a, const CInteger& b);
		friend CInteger operator % (const CInteger& a, const CInteger& b);
		friend bool operator == (const CInteger& a, const CInteger& b);

		friend bool operator < (const CInteger& a, const CInteger& b);
		friend bool operator > (const CInteger& a, const CInteger& b);
		friend bool operator >= (const CInteger& a, const CInteger& b);

		friend CInteger operator / (const CInteger& a, const CInteger& b);

		friend CInteger operator - (const CInteger& a);

		friend CInteger& operator += (CInteger& a, const CInteger& b);
		friend CInteger& operator -= (CInteger& a, const CInteger& b);

		friend CInteger lcm(const CInteger& a, const CInteger& b);
		friend CInteger gcd(const CInteger& a, const CInteger& b);
		friend CExtendedGCDResults<CInteger> extendedGCD(const CInteger& a, const CInteger& b);
		friend CExtendedExtendedGCDResults<CInteger> extendedExtendedGCD(const CInteger& a, const CInteger& b);
	protected:
		mpz_t data;

		void initialize() { mpz_init(data); };

	public:
		typedef CRational<CInteger> RationalType;

		CInteger() { initialize(); }
		CInteger(const mpz_t& iData) { initialize(); mpz_set(data, iData); }
		CInteger(const mpz_class& iData) { initialize(); mpz_set(data, iData.get_mpz_t()); }

		~CInteger() { mpz_clear(data); };

		mpz_class getMPZClass() const { return mpz_class(data); }

		mpz_t& getGMPData() { return data; }

		signed long int toInt() const { return mpz_get_si(data); }

		CInteger(long int i) {
			initialize();

			*this = i;
		}

		CInteger(const CInteger& i) {
			initialize();

			*this = i;
		}

		/*CInteger& operator = (unsigned long int i) {
			mpz_set_ui(data, i);

			return *this;
		};*/

		CInteger& operator = (long int i) {
			mpz_set_si(data, i);

			return *this;
		};

		CInteger& operator = (const CInteger& i) {
			if (this != &i) mpz_set(data, i.data);

			return *this;
		};

		void abs() { mpz_abs(data, data); }

		CInteger getAbs() const {
			CInteger result = CInteger();

			result = *this;
			result.abs();

			return result;
		}

		CInteger getLCMCofactor(CInteger& other) const;

		string toString() const { return string(mpz_get_str(NULL, 10, data)); }

		void print() const;

		bool isZero() const;
		bool isOne() const;

		static CInteger getZero() { return CInteger(0); }
		static CInteger getOne() { return CInteger(1); }

		string toStringHtml() const { return toString(); };

		static std::string getTypeName() { return "I"; }
		static std::string getTypeNameHtml() { return "I"; }
	};

	CInteger gcd(const CInteger& a, const CInteger& b);
	CExtendedGCDResults<CInteger> extendedGCD(const CInteger& a, const CInteger& b);
	CExtendedExtendedGCDResults<CInteger> extendedExtendedGCD(const CInteger& a, const CInteger& b);

	CInteger lcm(const CInteger& a, const CInteger& b);

	bool operator == (const CInteger& a, const CInteger& b);
	bool operator != (const CInteger& a, const CInteger& b);

	bool operator >= (const CInteger& a, const CInteger& b);

	CInteger operator + (const CInteger& a, const CInteger& b);
	CInteger operator - (const CInteger& a, const CInteger& b);
	CInteger operator * (const CInteger& a, const CInteger& b);
	CInteger operator % (const CInteger& a, const CInteger& b);

	bool operator < (const CInteger& a, const CInteger& b);
	bool operator > (const CInteger& a, const CInteger& b);

	bool sameSign(const CInteger& a, const CInteger& b);

	CInteger operator / (const CInteger& a, const CInteger& b);

	CInteger operator - (const CInteger& a);

	CInteger& operator += (CInteger& a, const CInteger& b);
	CInteger& operator -= (CInteger& a, const CInteger& b);
}

// #include "Rational.h"

#endif
