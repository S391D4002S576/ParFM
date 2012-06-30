#ifndef TEST_H_
#define TEST_H_

#include "../cute/cute.h"
#include "tests/Test_MultiListener.h"

#include "tests/Test_Funcs.h"

using cute::suite;

cute::Test_MultiListener test_MultiListener;

suite* Test_Integer_runSuite();
suite* Test_Vector_runSuite();
suite* Test_Matrix_runSuite();

suite* Test_ConnectedComponentFinder_runSuite();

suite* Test_Module_runSuite();
//suite* Test_Lattice_runSuite();
suite* Test_Lattice_NewHadron_runSuite();
suite* Test_Cone_runSuite();

suite* Test_ModuleRelation_runSuite();
suite* Test_LatticeRelation_runSuite();

suite* Test_FourierMotzkin_runSuite();
suite* Test_LimsLinearisation_runSuite();

suite* Test_IntegerProgramming_runSuite();

suite* Test_Graph_runSuite();
suite* Test_RADGraph_runSuite();
suite* Test_RDMGraph_runSuite();

suite* Test_DeskriptorSet_runSuite();

suite* Test_Computation_runSuite();
suite* Test_SetRelation_NewHadron_runSuite();
suite* Test_Computation_NewHadron_runSuite();

suite* Test_FinitePartition_runSuite();

#endif /* TEST_H_ */
