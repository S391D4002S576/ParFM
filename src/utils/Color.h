#ifndef COLOR_H_
#define COLOR_H_

#include <fstream>
using std::ostream;

class CColor {
		friend ostream& operator << (ostream &os, const CColor& color);
	public:
		int rgbValue;
		
		CColor(int iRGBValue) : rgbValue(iRGBValue) { };
		
		CColor& operator = (const CColor& other) {
			if (this != &other) {
				rgbValue = other.rgbValue;
			}
			
			return *this;
		}
};

ostream& operator << (ostream& os, const CColor& color);

#endif /* COLOR_H_ */
