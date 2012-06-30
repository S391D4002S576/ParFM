#ifndef SVGECIRCLE_H_
#define SVGECIRCLE_H_

#include "SVGElement.h"
#include "Color.h"
#include "Point.h"

class CSVGECircle : public CSVGElement {
	public:
		CPoint center;
		int radius;
		CColor fillColor;
		CColor strokeColor;
		int strokeWidth;
		
		CSVGECircle(CPoint iCenter, int iRadius, CColor iFillColor, CColor iStrokeColor, int iStrokeWidth)
		: center(iCenter), radius(iRadius), fillColor(iFillColor), strokeColor(iStrokeColor), strokeWidth(iStrokeWidth)
		{ };		

		void encode(ostream& os) {
			os 
			<< "<circle cx=\"" << center.x 
			<< "\" cy=\"" << center.y 
			<< "\" r=\"" << radius 
			<< "\" fill=\"" << fillColor
			<< "\" stroke=\"" << strokeColor
			<< "\" stroke-width=\"" << strokeWidth
			<< "\" />\n";
		}
		
		virtual ~CSVGECircle() { }
};

#endif /* SVGECIRCLE_H_ */
