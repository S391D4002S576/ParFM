#ifndef TOSTRING_CPP_
#define TOSTRING_CPP_

#include "ToString.h"

string intToStr(int x) {
	stringstream ss;
	ss << x;
	return ss.str();
}

#endif /* TOSTRING_CPP_ */
