#ifndef SVGELEMENT_CPP_
#define SVGELEMENT_CPP_

#include "SVGElement.h"

ostream& operator << (ostream& os, CSVGElement* svgElement) {
	svgElement->encode(os);
	
	return os;
}

#endif /* SVGELEMENT_CPP_ */
