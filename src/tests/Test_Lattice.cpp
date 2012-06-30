#ifdef OBSOLETE

#include "../cute/cute.h"
#include "../cute/cute_suite.h"

#include "../partyql/graph/Graph.h"
#include "../partyql/basicmath/scalar/Integer.h"
#include "../partyql/basicmath/Vektor.h"
#include "../partyql/basicmath/Lattice.h"
#include "../partyql/basicmath/Hadron.h"
#include "../partyql/basicmath/VectorDescriptor.h"

#include "../utils/SVG.h"
#include "../utils/SVGELattice.h"

#include "Test_Funcs.h"

using namespace AlgoTrans;

namespace Test_Lattice {


	/*void CLattice_AffineExemplar() {
		test_Lattice_initializeIntegerVectors();

		CLattice<I> a = CLattice<I>(2);
		a.addGeneratingVertex(new CVector<I>(iv_5_1));
		a.addGeneratingVertex(new CVector<I>(iv1_1));

		IVector* exemplar = a.getAffineExemplar();

		exemplar->print();

		ASSERT((exemplar != NULL) && gcd((*exemplar)[0] - 1, 6) == 6);

		delete exemplar;
	}

	void CLattice_AffineExemplar_NonUnitHomogenizedCoordInGenVectors() {
		test_Lattice_initializeIntegerVectors();

		CLattice<I> a = CLattice<I>(3);
		a.addGeneratingVertex(new CVector<I>(iv1_0_2));
		a.addGeneratingVertex(new CVector<I>(iv0_1_3));
		a.addGeneratingVertex(new CVector<I>(iv0_0_6));

		IVector* exemplar = a.getAffineExemplar();

		exemplar->print();

		ASSERT((exemplar != NULL) && gcd((*exemplar)[0] - 5, 6) == 6);
		ASSERT((exemplar != NULL) && gcd((*exemplar)[1] - 5, 6) == 6);

		delete exemplar;
	}*/

	void CLattice_AffineElementOf_1D() {
		CLattice<I> a = CLattice<I>(G, matrixToDescriptorSet(IMx(2,
											IV_(2, 1, -5), IV_(2, 1, 1))));

		ASSERT(a.containsAffinePoint(IV(1, -5)));
		ASSERT(a.containsAffinePoint(IV(1, 1)));
		ASSERT(a.containsAffinePoint(IV(1, 7)));
		ASSERT(a.containsAffinePoint(IV(1, 13)));

		ASSERT(!a.containsAffinePoint(IV(1, -4)));
		ASSERT(!a.containsAffinePoint(IV(1, -3)));
		ASSERT(!a.containsAffinePoint(IV(1, -2)));
		ASSERT(!a.containsAffinePoint(IV(1, -1)));
		ASSERT(!a.containsAffinePoint(IV(1, 0)));
		ASSERT(!a.containsAffinePoint(IV(1, 2)));
		ASSERT(!a.containsAffinePoint(IV(1, 3)));
		ASSERT(!a.containsAffinePoint(IV(1, 4)));
		ASSERT(!a.containsAffinePoint(IV(1, 5)));
		ASSERT(!a.containsAffinePoint(IV(1, 6)));
		ASSERT(!a.containsAffinePoint(IV(1, 8)));
	}

	void CLattice_AffineElementOf_2D() {
		CLattice<I> a = CLattice<I>(G, matrixToDescriptorSet(IMx(3,
		IV_(3, 0, 6, 1),
		IV_(3, 1, 5, 0),
		IV_(3, 0, 5, 2))));

		ASSERT(a.containsAffinePoint(IV(2, 5, 0)));
		ASSERT(a.containsAffinePoint(IV(2, 11, 1)));
		ASSERT(a.containsAffinePoint(IV(2, 6, -1)));

		ASSERT(!a.containsAffinePoint(IV(2, 4, 0)));
		ASSERT(!a.containsAffinePoint(IV(2, 12, 1)));
		ASSERT(!a.containsAffinePoint(IV(2, 7, -1)));
		ASSERT(!a.containsAffinePoint(IV(2, 5, 1)));
		ASSERT(!a.containsAffinePoint(IV(2, 11, 2)));
		ASSERT(!a.containsAffinePoint(IV(2, 6, -2)));

		CSVGELattice<I> svgLat = CSVGELattice<I>(a, IntV(2, 0, 0), IntV(2, 20, 20), 50000, true);
		CSVG(&svgLat, true).saveToFile("CLattice_AffineElementOf_2D.svg");
	}
}

cute::suite* Test_Lattice_runSuite(){
	cute::suite& s = *(new cute::suite("Lattice"));

	s.push_back(CUTE(Test_Lattice::CLattice_AffineElementOf_1D));
	s.push_back(CUTE(Test_Lattice::CLattice_AffineElementOf_2D));

	return &s;
}

#endif
