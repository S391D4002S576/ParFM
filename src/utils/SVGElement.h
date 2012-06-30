#ifndef SVGELEMENT_H_
#define SVGELEMENT_H_

#include <fstream>
using std::ostream;

class CSVGElement {
	friend ostream& operator << (ostream &os, CSVGElement* svgElement);
	
	public:
		virtual void encode(ostream& os) = 0;
		
		virtual ~CSVGElement() { };
};

ostream& operator << (ostream& os, CSVGElement* svgElement);

#endif /*SVGELEMENT_H_*/
