
#include "../cute/cute.h"
#include "../cute/ide_listener.h"

#include "../partyql/graph/Graph.h"
#include "../partyql/basicmath/scalar/Integer.h"

using namespace AlgoTrans;

namespace Test_Integer {
CInteger i2, i3, i5, i7, i10, i11, i21;
CInteger i2_, i3_, i5_, i7_, i10_, i11_, i21_;

void initializeIntegers() {
	i2 = CInteger(2);
	i3 = CInteger(3);
	i5 = CInteger(5);
	i7 = CInteger(7);
	i10 = CInteger(10);
	i11 = CInteger(11);
	i21 = CInteger(21);

	i2_ = CInteger(2);
	i3_ = CInteger(3);
	i5_ = CInteger(5);
	i7_ = CInteger(7);
	i10_ = CInteger(10);
	i11_ = CInteger(11);
	i21_ = CInteger(21);
}

void CInteger_equality() {
	initializeIntegers();

	ASSERT(i2 == i2_);
	ASSERT(i3 == i3_);
	ASSERT(i5 == i5_);
	ASSERT(i7 == i7_);
	ASSERT(i10 == i10_);
	ASSERT(i11 == i11_);
	ASSERT(i21 == i21_);

	ASSERT(!(i2 == i3_));
	ASSERT(!(i7 == i5_));
	ASSERT(!(i11 == i10));
	ASSERT(!(i21 == i3));
}

void CInteger_copy() {
	initializeIntegers();

	CInteger i = i2;
	ASSERT(i == i2_);
	ASSERT(i2 == i);
}

void CInteger_inequality() {
	ASSERT(i2 != i7_);
	ASSERT(i2 != i11_);
	ASSERT(i5 != i10);
	ASSERT(i21 != i11);

	ASSERT(!(i2 != i2_));
	ASSERT(!(i3 != i3_));
	ASSERT(!(i5 != i5_));
	ASSERT(!(i7 != i7_));
	ASSERT(!(i10 != i10_));
	ASSERT(!(i11 != i11_));
	ASSERT(!(i21 != i21_));
}

void CInteger_addition() {
	initializeIntegers();

	ASSERT(i2 + i3 == i5);
	ASSERT(i5 + i2 == i7);
	ASSERT(i7 + i3 == i10);
	ASSERT(i10 + i11 == i21);

	ASSERT(i7 + i7 + i7 == i21);
	ASSERT(i7 + i3 + i11 == i21);
}

void CInteger_multiplication() {
	initializeIntegers();

	ASSERT(i3 * i7 == i21);
	ASSERT(i5 * i2 == i10);
}

}

cute::suite* Test_Integer_runSuite(){
	cute::suite& s = *(new cute::suite("Integer"));

	Test_Integer::initializeIntegers();

	s.push_back(CUTE(Test_Integer::CInteger_equality));
	s.push_back(CUTE(Test_Integer::CInteger_copy));
	s.push_back(CUTE(Test_Integer::CInteger_inequality));
	s.push_back(CUTE(Test_Integer::CInteger_addition));
	s.push_back(CUTE(Test_Integer::CInteger_multiplication));

	return &s;
}


