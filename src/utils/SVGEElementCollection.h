#ifndef SVGEELEMENTCOLLECTION_H_
#define SVGEELEMENTCOLLECTION_H_

#include "SVGElement.h"

class CSVGEElementCollection  : public CSVGElement {
	public:
		vector<CSVGElement* > elements;
		
		void encode(ostream& os) {
			for (unsigned int q = 0; q < elements.size(); q++) {
				os << elements[q];
			}
		}
		
		virtual ~CSVGEElementCollection() {
			for (int q = elements.size() - 1; q >= 0; q--) {
				delete elements[q];
			}
		}
};

#endif /*SVGEELEMENTCOLLECTION_H_*/
