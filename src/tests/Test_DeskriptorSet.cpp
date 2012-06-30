#include "../cute/cute.h"

#include "Test_Funcs.h"

#include <fstream>
using std::ofstream;

#include "../partyql/core/scalar/Integer.h"
#include "../partyql/core/scalar/Bool.h"
#include "../partyql/core/Descriptor.h"

using namespace AlgoTrans;

namespace Test_DeskriptorSet {

	void hadron_CombineDeskriptorSet() {
		IHadron a = IHC(false, 1, IV_(2, 0, 1));
		IHadron b = IHC(false, 1, IV_(2, 1, 0));

		IHadron aANDb = IHC(false, 2, IV_(2, 1, 0), IV_(2, 0, 1));
		IHadron a_b = combineDeskriptorSets(a, b, C);

		// XXX: ASSERT(aANDb == a_b);
	}

	void flat_containment() {
		Polyheder hA = IP(C, 2, IV_(2, 1, 0), IV_(2, 0, 1));
		Polyheder hB = IP(G, 4, IV_(2, 1, 0), IV_(2, 0, 1), IV_(2, -1, 0), IV_(2, 0, -1));

		ASSERT(!(hB <= hA));
	}

	void flat_Equality() {
		Polyheder hA = IP(C, 2, IV_(2, 1, 0), IV_(2, 0, 1));
		Polyheder hB = IP(G, 4, IV_(2, 1, 0), IV_(2, 0, 1), IV_(2, -1, 0), IV_(2, 0, -1));

		ASSERT(hA != hB);
	}

	void deskriptorSet_nonNegativelyGenerates() {
		IDeskriptorSet a = IDS(false, 4, IV_(5, 0, -1, -3, -1,  1),
				                         IV_(5, 0,  1, -1, -3, -1),
				                         IV_(5, 0,  0,  0,  1,  1),
				                         IV_(5, 0,  1,  1,  1,  1));

        IDeskriptor d = IDeskriptor(IV(5, 0, 1, -1, 1,  3), false);

        if (!a.nonNegativelyGenerates(d)) {
        	Bool(a.nonNegativelyGenerates(d)).print();
        }
        ASSERT(a.nonNegativelyGenerates(d));
	}

	void deskriptorSet_nonNegativelyGenerates_2() {
		IDeskriptorSet a = IDS(false, 7, IV_(5, 0, -1, -3, -1,  1),
				                         IV_(5, 0,  1, -1, -3, -1),
				                         IV_(5, 0,  0,  0,  1,  1),
				                         IV_(5, 0,  1,  1,  1,  1),
				                         IV_(5, 2, -1,  1,  5,  3),
				                         IV_(5, 2,  1, -1, -1,  1),
				                         IV_(5, 2,  3,  5,  1, -1));

        IDeskriptor d = IDeskriptor(IV(5, 0, 1, -1, 1,  3), false);

        ASSERT(a.nonNegativelyGenerates(d));

	}

	void deskriptorSet_flatConcatenation() {
		/*IDeskriptorSet a = IDS(true, IV_(3, 1, 1, 0), IV_(3, 0, 0, 1));
		IDeskriptorSet b = IDS(true, IV_(2, -1, 0), IV_(2, 0, -1));

		 1,  1, 0, 0
		 0,  0, 1, 0
		-1,  0, 0, 0
		 0, -1, 0, 0
		C: DS<D<I>>(4, D<I>(  1 -1  0  0 == ))
		H<DS<D<I>>>[G: DS<D<I>>(4, D<I>( 0 0 0 1 == ), D<I>( 0 0 1 0 == ), D<I>( 1 1 0 0 == )) | C: DS<D<I>>(4, D<I>(  1 -1  0  0 == ))]
		*/
	}

	void deskriptorSet_Simplify() {
		IDeskriptorSet a = IDS(false, 3, IV_(2, 1, 0), IV_(2, 0, 1), IV_(2, -1, -1));
		IDeskriptorSet b = a.getSimplifiedDeskriptorSet();

		//b.print();

		ASSERT(a == b);
	}

	void deskriptorSet_Simplify2() {
		IDeskriptorSet a = IDS(false, 4, IV_(2, 1, 0), IV_(2, 0, 1), IV_(2, 1, 0), IV_(2, 0, 1));
		IDeskriptorSet b = a.getSimplifiedDeskriptorSet();
		IDeskriptorSet bExp = IDS(false, 2, IV_(2, 1, 0), IV_(2, 0, 1));

		//b.print();

		ASSERT(b == bExp);
		ASSERT(a == b);
	}

	/*void deskriptorSet_Simplify3() {
		IDeskriptorSet a = IDS(false, 4, IV_(6, 0, 0, 0, 0, 0, 1), IV_(6, 0, 0, 0, 0, 1, 0), IV_(6, 0, 0, 0, 1, 0, 0), IV_(6, 0, 0, 1, 0, 0, 0));
		IDeskriptorSet b = IDS(true, 2, IV_(6, 0, 1, 0, -1, 1, -1), IV_(6, 1, 0, -1, 0, -1, 1));
		IDeskriptorSet ab = a >> b;

		IDeskriptorSet result = ab.getSimplifiedDeskriptorSet();
		result.print();

		IDeskriptorSet aX = IDS(false, 4, IV_(5, 0, 0, 0, 0, 1), IV_(5, 0, 0, 0, 1, 0), IV_(5, 0, 0, 1, 0, 0), IV_(5,  1,  0,  0, -1,  1));
		IDeskriptorSet bX = IDS(true, 1, IV_(5, 0, 1, -1, 1, -1));
		IDeskriptorSet abX = aX >> bX;

		ASSERT(result == abX);
	}*/
	/*void deskriptor_Projection() {
		IDeskriptorSet a = IDS(true, 2, IV_(4, 1, 0, 1, -1), IV_(4, 0, 1, -1, 1));
		IDeskriptor b = IDeskriptor(IV(4, 0, 0, 0, 1), false);

		IDeskriptor c = b.getProjection(a);

		IDeskriptor cExp = IDeskriptor(IV(4, 0, 0, 0, 1), false);

		ASSERT(c == cExp);
	}*/

	void deskriptorSet_FMProjection() {
		IDeskriptorSet a = IDS(false, 3, IV_(2, 1, 0), IV_(2, 0, 1), IV_(2, -1, -1));
		IDeskriptorSet b = a.getFMProjection(0, true);
		IDeskriptorSet bExp = IDS(true, 1, IV_(1, 1));

		//b.print();

		ASSERT(b == bExp);
	}

	void deskriptorSet_FMProjection_2() {
		IDeskriptorSet a = IDS(false, 4, IV_(4, 0, 0, 0, 1), IV_(4, 0, 0, 1, 0), IV_(4, 0, 1, 0, 0), IV_(4, 1, 0, 0, 0));
		IDeskriptorSet b = IDS(true, 2, IV_(4, 0, -1, 1, -1), IV_(4, -1, 0, -1, 1));
		IDeskriptorSet ab = a >> b;

		//ab.getSimplifiedDeskriptorSet().print();

		IDeskriptorSet result = ab.getFMProjection(0, true);
		//result.print();

		IDeskriptorSet aX = IDS(false, 4, IV_(3, 0, 0, 1), IV_(3, 0, 1, 0), IV_(3, 1, 0, 0), IV_(3,  0, -1,  1));
		IDeskriptorSet bX = IDS(true, 1, IV_(3, -1, 1, -1));
		IDeskriptorSet abX = aX >> bX;

		//abX.getSimplifiedDeskriptorSet().print();

		//ASSERT(result == abX);
	}

	void deskriptorSet_FMProjection_3() {
		IDeskriptorSet a = IDS(false, 4, IV_(6, 0, 0, 0, 0, 0, 1), IV_(6, 0, 0, 0, 0, 1, 0), IV_(6, 0, 0, 0, 1, 0, 0), IV_(6, 0, 0, 1, 0, 0, 0));
		IDeskriptorSet b = IDS(true, 2, IV_(6, 0, 1, 0, -1, 1, -1), IV_(6, 1, 0, -1, 0, -1, 1));
		IDeskriptorSet ab = a >> b;

		//ab.getSimplifiedDeskriptorSet().print("ab");

		IDeskriptorSet result = ab.getFMProjection(2, true);
		//result.print("result");

		IDeskriptorSet aX = IDS(false, 4, IV_(5, 0, 0, 0, 0, 1), IV_(5, 0, 0, 0, 1, 0), IV_(5, 0, 0, 1, 0, 0), IV_(5,  1,  0,  0, -1,  1));
		IDeskriptorSet bX = IDS(true, 1, IV_(5, 0, 1, -1, 1, -1));
		IDeskriptorSet abX = aX >> bX;

		ASSERT(result == abX);
	}

	void deskriptorSet_FMProjection_4() {
		IDeskriptorSet a = IDS(true, 2, IV_(4, 0, 1, 0, -1), IV_(4, 1, 0, 0, -1));
		IDeskriptorSet aX = IDS(true, 1, IV_(3, 1, -1, 0));

		IDeskriptorSet result = a.getFMProjection(3, true);
		//result.print("result");

		ASSERT(result == aX);
	}

	void deskriptorSet_Dual_A() {
		IDeskriptorSet a = IDS(false, 4, IV_(2, 1, 0), IV_(2, 0, 1), IV_(2, 1, -1), IV_(2, -1, 1));
		IDeskriptorSet aDual = a.getDual();
		//aDual.print("aDual");

		IDeskriptorSet b = IDS(false, 1, IV_(2, 1, 1));

		ASSERT(aDual == b);
	}

	void deskriptorSet_Dual_B() {
		IDeskriptorSet a = IDS(true, 2, IV_(3, 1, 1, 0), IV_(3, 0, 0, 1));
		IDeskriptorSet aDual = a.getDual();
		//aDual.print("aDual");

		IDeskriptorSet b = IDS(true, 1, IV_(3, 1, -1, 0));

		ASSERT(aDual == b);
	}

	void deskriptorSet_Dual() {
		IDeskriptorSet a = IDS(false, 5, IV_(3, 1, 0, 0), IV_(3, 0, 1, 0), IV_(3, 0, 0, 1), IV_(3, 0, 1, -1), IV_(3, 0, -1, 1));
		IDeskriptorSet aDual = a.getDual();
		//aDual.print();

		IDeskriptorSet b = IDS(false, 2, IV_(3, 1, 0, 0), IV_(3, 0, 1, 1));

		ASSERT(aDual == b);
	}

	void polyheder_Projection() {
		Polyheder a = IP(C, 3, IV_(2, 1, 0), IV_(2, 0, 1), IV_(2, -1, -1));
		Polyheder aProj = a.getProjection(std::vector<int>(1, 0));

		//aProj.print();
		Polyheder b = IP(C, 2, IV_(1, 1), IV_(1, -1));

		ASSERT(aProj == b);
	}

	void polyheder_Empty() {
		Polyheder a = IPC(false, 3, IV_(2, 1, 0), IV_(2, 0, 1), IV_(2, -1, -1));

		ASSERT(a.isEmpty());
	}

	void polyheder_Join() {
		Polyheder a = IPC(false, 1, IV_(2, 0, 1));
		Polyheder b = IPC(false, 1, IV_(2, 1, 0));

		Polyheder aANDb = IPC(false, 2, IV_(2, 1, 0), IV_(2, 0, 1));
		Polyheder a_b = a && b;
		//aANDb.print();
		//a_b.print();

		ASSERT_EQ(aANDb, a_b);
	}

	void polyhedralDomain_Join() {
		PolyhedralDomain a = IPDC(false, 1, IV_(2, 0, 1));
		PolyhedralDomain b = IPDC(false, 1, IV_(2, 1, 0));

		PolyhedralDomain aANDb = IPDC(false, 2, IV_(2, 1, 0), IV_(2, 0, 1));
		PolyhedralDomain a_b = a && b;

		ASSERT_EQ(aANDb, a_b);
	}

	void polyheder_Empty2() {
		Polyheder a = IPC(false, 1, IV_(3,  1,  1, -1)) && IPC(false, 2, IV_(3,  0,  1,  0), IV_(3,  1,  0,  0));

		ASSERT(!a.isEmpty());
	}

	void deskriptorSet_Equality() {
		IDeskriptorSet a = IDS(true, 3, IV_(4, 0, 0, 1, 1), IV_(4, 0, 1, 0, 0), IV_(4, 1, 0, 0, 0)) + IDS(false, 1, IV_(4, 0, 0, 0, 1));
		IDeskriptorSet b = IDS(true, 3, IV_(4, 0, 0, 1, 1), IV_(4, 0, 1, 0, 0), IV_(4, 1, 0, 0, 0)) + IDS(false, 1, IV_(4, 1, 1, 0, 1));

		ASSERT_EQ(a, b);
	}

	void dimensionality1() {
		IHadron h = IHadron(C, IDS(true, 1, IV_(3, 0, 1, -1)) + IDS(false, 2, IV_(3, 0, 1, 0), IV_(3, 1, 0, 0)));
		CInteger(h.getDimensionality()).print();
		CInteger(h.getDimensionality()).print();
		ASSERT(h.getDimensionality() == 2);
	}
}

cute::suite* Test_DeskriptorSet_runSuite(){
	cute::suite& s = *(new cute::suite("DeskriptorSet"));

	s.push_back(CUTE(Test_DeskriptorSet::deskriptorSet_nonNegativelyGenerates));
	s.push_back(CUTE(Test_DeskriptorSet::deskriptorSet_nonNegativelyGenerates_2));

//	s.push_back(CUTE(Test_DeskriptorSet::deskriptor_Projection));

	s.push_back(CUTE(Test_DeskriptorSet::deskriptorSet_Simplify));
	s.push_back(CUTE(Test_DeskriptorSet::deskriptorSet_Simplify2));
//	s.push_back(CUTE(Test_DeskriptorSet::deskriptorSet_Simplify3));

	s.push_back(CUTE(Test_DeskriptorSet::deskriptorSet_FMProjection));
	s.push_back(CUTE(Test_DeskriptorSet::deskriptorSet_FMProjection_2));
	s.push_back(CUTE(Test_DeskriptorSet::deskriptorSet_FMProjection_3));
	s.push_back(CUTE(Test_DeskriptorSet::deskriptorSet_FMProjection_4));

	s.push_back(CUTE(Test_DeskriptorSet::deskriptorSet_Dual_A));
	s.push_back(CUTE(Test_DeskriptorSet::deskriptorSet_Dual));
	s.push_back(CUTE(Test_DeskriptorSet::deskriptorSet_Dual_B));

	s.push_back(CUTE(Test_DeskriptorSet::deskriptorSet_Equality));

	s.push_back(CUTE(Test_DeskriptorSet::hadron_CombineDeskriptorSet));

	s.push_back(CUTE(Test_DeskriptorSet::flat_containment));
	s.push_back(CUTE(Test_DeskriptorSet::flat_Equality));

	s.push_back(CUTE(Test_DeskriptorSet::polyheder_Projection));
	s.push_back(CUTE(Test_DeskriptorSet::polyheder_Empty));
	s.push_back(CUTE(Test_DeskriptorSet::polyheder_Empty2));
	s.push_back(CUTE(Test_DeskriptorSet::polyheder_Join));
	s.push_back(CUTE(Test_DeskriptorSet::polyhedralDomain_Join));

	s.push_back(CUTE(Test_DeskriptorSet::dimensionality1));


	return &s;
}

