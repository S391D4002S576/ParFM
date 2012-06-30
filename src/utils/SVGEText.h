#ifndef SVGETEXT_H_
#define SVGETEXT_H_

#include <string>
using std::string;

#include "SVGElement.h"
#include "Color.h"
#include "Point.h"

class CSVGEText : public CSVGElement {
	public:
		CPoint location;
		CColor color;
		string text;
		int fontSize;
		string fontFamily;
		string textAnchor;
		string dominantBaseline;
		bool vertical;
		
		CSVGEText(CPoint iLocation, CColor iColor, string iText, int iFontSize, string iFontFamily, string iTextAnchor, string iDominantBaseline)
		: location(iLocation), color(iColor), text(iText), fontSize(iFontSize), fontFamily(iFontFamily), textAnchor(iTextAnchor), dominantBaseline(iDominantBaseline)
		{ vertical = false; };		

		void encode(ostream& os) {
			os 
			<< "<text x=\"" << location.x 
			<< "\" y=\"" << location.y 
			<< "\" font-family=\"" << fontFamily 
			<< "\" font-size=\"" << fontSize 
			<< "\" fill=\"" << color
			<< "\" text-anchor=\"" << textAnchor
			<< "\" dominant-baseline=\"" << dominantBaseline << "\""
			;
			
			if (vertical) {
				os << " transform=\"rotate(90, " << location.x << ", " << location.y << ")\"";
			}
			
			os
			<< ">"
			<< text
			<< "</text>\n";
		}
		
		virtual ~CSVGEText() { }
};

#endif /* SVGECIRCLE_H_ */
