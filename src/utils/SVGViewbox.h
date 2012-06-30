#ifndef SSSVGVIEWBOX_H_
#define SSSVGVIEWBOX_H_

#include <vector>
using std::vector;
using std::ios_base;

#include <iostream>
#include <fstream>
using std::ostream;
using std::ofstream;

#include <string>
using std::string;

#include "Point.h"
#include "Color.h"

#include "SVGElement.h"
#include "SVGEElementCollection.h"

#include "SVGELine.h" 
#include "SVGECircle.h"

class CSVGEViewBox: public CSVGElement {
	public:
		CSVGElement* subElement;
		int xPct;
		int yPct;
		int wPct;
		int hPct;
		int width;
		int height;
		
		CSVGEViewBox(CSVGElement* iSubElement, int iXPct, int iYPct, int iWPct, int iHPct, int iWidth, int iHeight) 
		: subElement(iSubElement), xPct(iXPct), yPct(iYPct), wPct(iWPct), hPct(iHPct), width(iWidth), height(iHeight) {
		}
		
		void encode(ostream& os) {
			os
			<< "<svg x=\"" << xPct 
			<< "\" y=\"" << yPct 
			<< "\" width=\"" << wPct 
			<< "\" height=\"" << hPct  
			<< "\" viewBox=\"0 0 " << width << " " << height << "\">\n"
			<< subElement
			<< "</svg>";	
		}
		
		virtual ~CSVGEViewBox() { delete subElement; };
};

#endif /*SSSVGVIEWBOX_H_*/
