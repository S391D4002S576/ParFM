#include "../cute/cute.h"
#include "../cute/ide_listener.h"
#include "../cute/cute_suite.h"

#include "Test_Funcs.h"

#include "../partyql/graph/Graph.h"
#include "../partyql/basicmath/scalar/Integer.h"
#include "../partyql/basicmath/Vektor.h"
#include "../partyql/basicmath/Module.h"

using namespace AlgoTrans;

namespace Test_Module {
	typedef CInteger I;

	// Depends on HNF reduction through Module-equality
	void CModule_Orthogonalisation() {
		// TEST: We might additionally want to test that the resulting vertices are "unimodular" (e.g. gcd = 1)
		bool verbose = false;

		// This one should have the same solution
		CModule<I> a = CModule<I>(G, matrixToDescriptorSet(IMx(3, IV_(3, 6, 0, 1),
											IV_(3, 2, -3, -1), IV_(3, 4, 3, 2))));
		a = !a;

		if (verbose) printf("a:\n%s\n", a.toString().c_str());

		CModule<I> exp = CModule<I>(G, matrixToDescriptorSet(IMx(1, IV_(3, 3, 8, -18))));

		if (verbose) printf("exp:\n%s", exp.toString().c_str());

		ASSERT(a == exp);

		// This one should have the same solution
		CModule<I> c = CModule<I>(G, matrixToDescriptorSet(IMx(3, IV_(3, 6, 0, 1),
											IV_(3, 6, 0, 1), IV_(3, 4, 3, 2))));
		c = !c;

		if (verbose) printf("c:\n%s\n", c.toString().c_str());

		ASSERT(c == exp);


		// Empty solution for this one
		CModule<I> fexp = CModule<I>(3, G);

		if (verbose) printf("exp:\n%s", fexp.toString().c_str());

		CModule<I> d = CModule<I>(G, matrixToDescriptorSet(IMx(3, IV_(3, 4, 4, 1),
				           IV_(3, 6, 0, 1), IV_(3, 4, 3, 2))));
		d = !d;

		if (verbose) printf("d:\n%s\n", d.toString().c_str());

		ASSERT(d == fexp);

		// This one should have the same solution
		CModule<I> e = CModule<I>(G, matrixToDescriptorSet(IMx(3, IV_(3, 4, 4, 1), IV_(3, 2, -3, -1), IV_(3, 4, 3, 2))));
		e = !e;

		if (verbose) printf("e:\n%s\n", e.toString().c_str());

		ASSERT(e == fexp);
	}

	/* Depends on HNF reduction through Module-equality
	 * This test failed previously due to forgetting to swap rows for the U matrix
	 */
	void CModule_Orthogonalisation_2() {
		CModule<I> a = CModule<I>(G, matrixToDescriptorSet(IMx(2, IV_(4, 1, 0, 0, 1), IV_(4, 0, 1, -1, 0))));
		a = !a;

		CModule<I> exp = CModule<I>(G, matrixToDescriptorSet(IMx(2, IV_(4, 1, 0, 0, -1), IV_(4, 0, 1, 1, 0))));

		ASSERT(a == exp);
	}



	void CModule_Intersection() {
		CModule<I> a = CModule<I>(G, matrixToDescriptorSet(IMx(3, IV_(4, 1, 0, 0, 1), IV_(4, 0, 1, 0, 0), IV_(4, 0, 0, 1, 0))));
		CModule<I> b = CModule<I>(G, matrixToDescriptorSet(IMx(2, IV_(4, 1, 0, 0, 1), IV_(4, 0, 1, -1, 0))));
		CModule<I> c = CModule<I>(G, matrixToDescriptorSet(IMx(2, IV_(4, 1, 0, 0, 1), IV_(4, 0, 1, -1, 0))));

		ASSERT(c == (a * b));
		ASSERT(c == (b * a));
		ASSERT(c == ((b * a) * c));
		ASSERT(c == ((b * a) * a));
		ASSERT(c == ((b * a) * b));
		ASSERT(c == (a * c));
		ASSERT(c == (b * c));
		ASSERT(c == ((a * a) * (b * b)));
		ASSERT(c == (((a * a) * (b * b)) * c));
		ASSERT(c == ((a * a) * ((b * b) * c)));

		a = CModule<I>(G, matrixToDescriptorSet(IMx(3, IV_(4, 1, 0, 0, 1), IV_(4, 0, 1, 0, 0), IV_(4, 0, 0, 1, 0))));
		b = CModule<I>(G, matrixToDescriptorSet(IMx(2, IV_(4, 1, 0, 0, -1), IV_(4, 0, 1, -1, 0))));
		c = CModule<I>(G, matrixToDescriptorSet(IMx(1, IV_(4, 0, 1, -1, 0))));

		ASSERT(c == (a * b));
		ASSERT(c == (b * a));
		ASSERT(c == ((b * a) * c));
		ASSERT(c == ((b * a) * a));
		ASSERT(c == ((b * a) * b));
		ASSERT(c == (a * c));
		ASSERT(c == (b * c));
		ASSERT(c == ((a * a) * (b * b)));
		ASSERT(c == (((a * a) * (b * b)) * c));
		ASSERT(c == ((a * a) * ((b * b) * c)));

		a = CModule<I>(G, matrixToDescriptorSet(IMx(2, IV_(4, 1, 0, 0, 1), IV_(4, 0, 1, 1, 0))));
		b = CModule<I>(G, matrixToDescriptorSet(IMx(2, IV_(4, 1, 0, 0, -1), IV_(4, 0, 1, -1, 0))));
		c = CModule<I>(4, G);

		ASSERT(c == (a * b));
		ASSERT(c == (b * a));
		ASSERT(c == ((b * a) * c));
		ASSERT(c == ((b * a) * a));
		ASSERT(c == ((b * a) * b));
		ASSERT(c == (a * c));
		ASSERT(c == (b * c));
		ASSERT(c == ((a * a) * (b * b)));
		ASSERT(c == (((a * a) * (b * b)) * c));
		ASSERT(c == ((a * a) * ((b * b) * c)));
	}
}

cute::suite* Test_Module_runSuite(){
	cute::suite &s = *(new cute::suite("Module"));

	s.push_back(CUTE(Test_Module::CModule_Orthogonalisation));
	s.push_back(CUTE(Test_Module::CModule_Orthogonalisation_2));
	s.push_back(CUTE(Test_Module::CModule_Intersection));

	return &s;
}

