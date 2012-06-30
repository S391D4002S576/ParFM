#include "../cute/cute.h"
#include "../cute/ostream_listener.h"
#include "../cute/ide_listener.h"
#include "../cute/cute_counting_listener.h"
#include "../cute/cute_suite.h"

#include "Test.h"

#include <sys/time.h>

using namespace AlgoTrans;

int main() {
	cute::Test_MultiListener multiListener;
	multiListener.addListener(new cute::ostream_listener());
	//multiListener.addListener(new cute::ide_listener());
	multiListener.addListener(new cute::counting_listener());

	cute::suite allSuite("");

	allSuite.push_back(Test_Integer_runSuite());

	allSuite.push_back(Test_Vector_runSuite());
	allSuite.push_back(Test_Matrix_runSuite());

	allSuite.push_back(Test_ConnectedComponentFinder_runSuite());

	allSuite.push_back(Test_Module_runSuite());
	//allSuite.push_back(Test_Lattice_runSuite());
	allSuite.push_back(Test_Cone_runSuite());

	allSuite.push_back(Test_Lattice_NewHadron_runSuite());

	allSuite.push_back(Test_ModuleRelation_runSuite());
	allSuite.push_back(Test_LatticeRelation_runSuite());

	allSuite.push_back(Test_FourierMotzkin_runSuite());
	// XXX: Replace Polyheder with cone in lims lin
	//allSuite.push_back(Test_LimsLinearisation_runSuite());

	allSuite.push_back(Test_IntegerProgramming_runSuite());

	allSuite.push_back(Test_Graph_runSuite());
	//allSuite.push_back(Test_FinitePartition_runSuite());

	allSuite.push_back(Test_RADGraph_runSuite());
	//allSuite.push_back(Test_RDMGraph_runSuite());

	allSuite.push_back(Test_DeskriptorSet_runSuite());

	allSuite.push_back(Test_SetRelation_NewHadron_runSuite());

	allSuite.push_back(Test_Computation_NewHadron_runSuite());
//	allSuite.push_back(Test_Computation_runSuite());

	int exitStatus = -1;

	double tPrev = rtclock();

	try { exitStatus = allSuite.performTest(&multiListener) ? 0 : 1; }
	catch (...) {
		printf("error: Test caused unexpected exception.\n");
	};

	double timeTests = (rtclock() - tPrev);
	printf("Test time = %0.6lfs\n",  timeTests);

	printTimeCollectors();

	return exitStatus;
}



