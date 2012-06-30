#ifndef ECHELONREDUCTION_H_
#define ECHELONREDUCTION_H_

#include <list>
using std::list;

#include "../scalar/ExtendedGCD.h"
#include "../Module.h"

namespace AlgoTrans {
	template <class R> class CModule;

	template <class R> class CEchelonReducer_VersionA {
	private:
		int rank;
		CModule<R>& X;
		list<int> pivotColumns;
	public:
		CEchelonReducer_VersionA(CModule<R>& iModule): X(iModule) { reduce(); }

		void reduce();
	};

	template <class R> class CEchelonReducer_VersionC {
		int rank;
		CModule<R>& X;
		list<int> pivotColumns;
	public:
		CEchelonReducer_VersionC(CModule<R>& iModule): X(iModule) { reduce(); }

		void reduce();
	};

	/* Does not give the same result as the other versions, but is that really a problem?
	 */
	template <class R> class CEchelonReducer_VersionD {
		int rank;
		CModule<R>& X;
		list<int> pivotColumns;
	public:
		CEchelonReducer_VersionD(CModule<R>& iModule): X(iModule) { reduce(); }

		void reduce();
	};

	template <class R> class CScaledEchelonReducer {
		int rank;
		CModule<R>& X;
		list<int> pivotColumns;
	public:
		CScaledEchelonReducer(CModule<R>& iModule): X(iModule) { reduce(); }

		void reduce();
	};
}

#endif /*ECHELONREDUCTION_H_*/
