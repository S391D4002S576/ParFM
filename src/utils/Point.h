#ifndef POINT_H_
#define POINT_H_

class CPoint {
	public:
		int x, y;
		
		CPoint(int iX, int iY) : x(iX), y(iY) { };
		
		CPoint& operator = (const CPoint& other) {
			if (this != &other) {
				x = other.x;
				y = other.y;
			}
			
			return *this;
		}
};

#endif /*POINT_H_*/
