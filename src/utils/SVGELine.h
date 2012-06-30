#ifndef SVGELINE_H_
#define SVGELINE_H_

#include "SVGElement.h"
#include "Color.h"
#include "Point.h"

class CSVGELine : public CSVGElement {
	public:
		CPoint from;
		CPoint to;
		CColor strokeColor;
		int strokeWidth;
		
		CSVGELine(CPoint iFrom, CPoint iTo, CColor iStrokeColor, int iStrokeWidth)
		: from(iFrom), to(iTo), strokeColor(iStrokeColor), strokeWidth(iStrokeWidth)
		{ };
		
		void encode(ostream& os) {
			os 
			<< "<line x1=\"" << from.x
			<< "\" y1=\"" << from.y
			<< "\" x2=\"" << to.x
			<< "\" y2=\"" << to.y
			<< "\" stroke=\"" << strokeColor
			<< "\" stroke-width=\"" << strokeWidth
			<< "\" />\n";
		}
		
		virtual ~CSVGELine() { };
};

#endif /*SVGELINE_H_*/
