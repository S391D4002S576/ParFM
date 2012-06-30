#include "Time.h"
#include "../../../../cute/cute.h"

double rtclock() {
    struct timezone Tzp;
    struct timeval Tp;
    int stat;
    stat = gettimeofday(&Tp, &Tzp);
    ASSERT(stat == 0); // error otherwise
    return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}

void printTimeCollectors() {
	printf("Times: -------------\n");
	for (unsigned int q = 0; q < timeCollectors.size(); q++) {
		TimeCollector& tc = *timeCollectors[q];
		printf("%s: %f s [%i%%] (%i) = %f s \n", tc.name.c_str(), tc.totalTime, (int) (100*(tc.totalTime/totalTimeCollected)), tc.runs, tc.totalTime/tc.runs);
	}
	printf("--------------------\n");
}

void resumeTimer(TimeCollector*& tc, std::string Name) {
	if (tc == NULL) {
		tc = new TimeCollector(Name);
	}
	tc->resume();
}

void pauseTimer(TimeCollector* tc) {
	tc->pause();
}

std::vector<TimeCollector*> timeCollectors;
double totalTimeCollected = 0.0;

TimeCollector * tcComputationClosureSpace;
TimeCollector * tcComputationClosureTime;

TimeCollector* tcDualisation;
TimeCollector* tcRDG;
TimeCollector* tcConstrainingProjection;
TimeCollector* tcEmpty;
TimeCollector* tcMeet;
TimeCollector* tcFastMeet;
TimeCollector* tcJoin;
TimeCollector* tcCRDC;
TimeCollector* tcCRUH;
TimeCollector* tcDDBidHull;
TimeCollector* tcDDSum;
TimeCollector* tcConcatOperation = NULL;
TimeCollector* tcFlatConcatOperation = NULL;
TimeCollector* tcFlatConcatOperationMeet = NULL;
TimeCollector* tcFlatConcatOperationProj = NULL;

TimeCollector* tcDomConeProjection = NULL;

TimeCollector* tcFullItSpace = NULL;
TimeCollector* tcUniHull = NULL;
TimeCollector* tcDeskriptors = NULL;

TimeCollector* tcComputationCalculateSpaceRDG = NULL;
TimeCollector* tcComputationCalculateTimeRDG = NULL;
TimeCollector* tcComputationCalcMultiTupleRelation = NULL;
TimeCollector* tcRDC2 = NULL;
TimeCollector* tcRDCBidHull = NULL;
TimeCollector* tcReferenceSpace = NULL;
TimeCollector* tcDomainCone = NULL;
TimeCollector* tcReflexiveDependenceCone = NULL;
TimeCollector* tcGetFullIterationSpace = NULL;


