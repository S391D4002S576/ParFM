#ifndef EXTENDEDGCD_H_
#define EXTENDEDGCD_H_

template <class R> class CExtendedExtendedGCDResults {
public:
	R leftBezoutFactor;
	R rightBezoutFactor;
	R leftLCMCofactor;
	R rightLCMCofactor;
	R gcd;
};

template <class R> class CExtendedGCDResults {
public:
	R leftBezoutFactor;
	R rightBezoutFactor;
	R gcd;
};

#endif /*EXTENDEDGCD_H_*/
