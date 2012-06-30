#ifndef COLOR_CPP_
#define COLOR_CPP_

#include "Color.h"

ostream& operator << (ostream& os, const CColor& color) {
	char hex[] = "0123456789ABCDEF";

	os << "#";
	for (int q = 5; q >= 0; q--) os << (hex[(color.rgbValue >> (4*q)) & 0xF]);
	
	return os;
}

#endif /* COLOR_CPP_ */
