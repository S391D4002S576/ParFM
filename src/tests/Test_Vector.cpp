#include "../cute/cute.h"
#include "../cute/ide_listener.h"
#include "../cute/cute_suite.h"

#include "../partyql/graph/Graph.h"
#include "../partyql/basicmath/scalar/Integer.h"
#include "../partyql/basicmath/Vektor.h"

using namespace AlgoTrans;

CInteger aiv3_5_7[] = {CInteger(3), CInteger(5), CInteger(7)};
CInteger aiv9_6_4[] = {CInteger(9), CInteger(6), CInteger(4)};
CInteger aiv12_11_11[] = {CInteger(12), CInteger(11), CInteger(11)};
CVector<CInteger> iv3_5_7, iv9_6_4, iv12_11_11;
CVector<CInteger> iv3_5_7_, iv9_6_4_, iv12_11_11_;

void initializeIntegerVectors() {
	iv3_5_7 = CVector<CInteger>(3, aiv3_5_7);
	iv3_5_7_ = CVector<CInteger>(3, aiv3_5_7);
	iv9_6_4 = CVector<CInteger>(3, aiv9_6_4);
	iv9_6_4_ = CVector<CInteger>(3, aiv9_6_4);
	iv12_11_11 = CVector<CInteger>(3, aiv12_11_11);
	iv12_11_11_ = CVector<CInteger>(3, aiv12_11_11);
}

void CVector_equality() {
	initializeIntegerVectors();

	ASSERT(iv3_5_7 == iv3_5_7);
	ASSERT(iv3_5_7 == iv3_5_7_);
	ASSERT(iv9_6_4 == iv9_6_4);
	ASSERT(iv9_6_4 == iv9_6_4_);

	ASSERT(!(iv3_5_7 == iv9_6_4_));
	ASSERT(!(iv9_6_4_ == iv12_11_11_));
	ASSERT(!(iv3_5_7 == iv12_11_11));
}

void CVector_copy() {
	initializeIntegerVectors();

	CVector<CInteger> v = iv3_5_7;

	ASSERT(v == iv3_5_7);
}

void CVector_inequality() {
	initializeIntegerVectors();

	ASSERT(iv3_5_7 != iv9_6_4_);
	ASSERT(iv9_6_4_ != iv12_11_11_);
	ASSERT(iv3_5_7 != iv12_11_11);
}

void CVector_addition() {
	initializeIntegerVectors();

	ASSERT(iv3_5_7 + iv9_6_4 == iv12_11_11);
	ASSERT(iv3_5_7 + iv9_6_4 == iv3_5_7_ + iv9_6_4_);
	ASSERT(iv3_5_7 + iv9_6_4 != iv12_11_11 + iv9_6_4_);
}

cute::suite* Test_Vector_runSuite(){
	cute::suite &s = *(new cute::suite("Vector"));

	initializeIntegerVectors();

	s.push_back(CUTE(CVector_equality));
	s.push_back(CUTE(CVector_copy));
	s.push_back(CUTE(CVector_inequality));
	s.push_back(CUTE(CVector_addition));
	//s.push_back(CUTE(CVector_multiplication));

	return &s;
}
