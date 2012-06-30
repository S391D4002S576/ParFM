#ifndef TIME_H_
#define TIME_H_

#include <sys/time.h>
#include <string>
#include <vector>
#include <stdio.h>

double rtclock();

class TimeCollector;

extern std::vector<TimeCollector*> timeCollectors;
extern double totalTimeCollected;

class TimeCollector {
public:
	double totalTime;
	int runs;
	std::string name;
	double lastTime;

	TimeCollector(std::string iName) : totalTime(0.0), runs(0), name(iName) { };

    void resume() { lastTime = rtclock(); runs++; };
    void pause() { double D = (rtclock() - lastTime); totalTime += D; totalTimeCollected += D; };
};

extern TimeCollector* tcDualisation;
extern TimeCollector* tcRDG;
extern TimeCollector* tcConstrainingProjection;
extern TimeCollector* tcEmpty;
extern TimeCollector* tcMeet;
extern TimeCollector* tcFastMeet;
extern TimeCollector* tcJoin;
extern TimeCollector* tcCRDC;
extern TimeCollector* tcCRUH;
extern TimeCollector* tcDDBidHull;
extern TimeCollector* tcDDSum;
extern TimeCollector* tcConcatOperation;
extern TimeCollector* tcFlatConcatOperation;
extern TimeCollector* tcFlatConcatOperationMeet;
extern TimeCollector* tcFlatConcatOperationProj;
extern TimeCollector * tcComputationClosureSpace;
extern TimeCollector * tcComputationClosureTime;

extern TimeCollector* tcDomConeProjection;

extern TimeCollector* tcFullItSpace;
extern TimeCollector* tcUniHull;
extern TimeCollector* tcDeskriptors;
extern TimeCollector* tcComputationCalculateSpaceRDG;
extern TimeCollector* tcComputationCalculateTimeRDG;
extern TimeCollector* tcComputationCalcMultiTupleRelation;
extern TimeCollector* tcRDC2;
extern TimeCollector* tcRDCBidHull;
extern TimeCollector* tcReferenceSpace;
extern TimeCollector* tcDomainCone;
extern TimeCollector* tcReflexiveDependenceCone;
extern TimeCollector* tcGetFullIterationSpace;



void resumeTimer(TimeCollector*& tc, std::string Name);
void pauseTimer(TimeCollector* tc);

void printTimeCollectors();

#define TIME(A, B) { \
	tPrev = rtclock(); \
	A \
	tDuration = (rtclock() - tPrev); tPrev = rtclock(); printf("T(%s) = %0.6lfs\n", B, tDuration); \
}

#define TIME_(A, tPrev, tDuration) { \
	tPrev = rtclock(); \
	A \
	tDuration = (rtclock() - tPrev); \
}

#endif
