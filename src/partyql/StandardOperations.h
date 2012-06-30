#ifndef STANDARDOPERATIONS_H_
#define STANDARDOPERATIONS_H_

namespace AlgoTrans {
	template <class X>
	class CShortestPathOperations {
		public:
			static X calcTransition(const X& a, const X& b) { return a + b; }

			static X calcCombination(const X& a, const X& b) { return (a < b) ? a : b; }
	};

	template <class X>
	class PathExistenceOperations {
		public:
			static X calcTransition(X a, X b) { return a && b; }

			static X calcCombination(X a, X b) { return a || b; }
	};
}

#endif /* STANDARDOPERATIONS_H_ */
