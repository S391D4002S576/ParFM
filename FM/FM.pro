LIBS += -lgmpxx -lgmp

HEADERS += \
    ../src/Test.h \
    ../src/partyql/Structures.h \
    ../src/partyql/Statement.h \
    ../src/partyql/StandardOperations.h \
    ../src/partyql/SpacePartitioner.h \
    ../src/partyql/RDMGraph.h \
    ../src/partyql/RADGraph.h \
    ../src/partyql/Polynomials.h \
    ../src/partyql/HierarPart.h \
    ../src/partyql/HierarGraph.h \
    ../src/partyql/Declarations.h \
    ../src/partyql/Computation.h \
    ../src/partyql/CommunicationChannel.h \
    ../src/partyql/basicmath/Vektor.h \
    ../src/partyql/basicmath/VectorDescriptor.h \
    ../src/partyql/basicmath/UnitaryTransformationDescriptor.h \
    ../src/partyql/basicmath/SetRelation.h \
    ../src/partyql/basicmath/ppl.h \
    ../src/partyql/basicmath/PolyhedralDomain.h \
    ../src/partyql/basicmath/Polyheder.h \
    ../src/partyql/basicmath/PartitionRelation.h \
    ../src/partyql/basicmath/Module.h \
    ../src/partyql/basicmath/Matrix.h \
    ../src/partyql/basicmath/Lattice.h \
    ../src/partyql/basicmath/IntegerProgramming.h \
    ../src/partyql/basicmath/HadronDescription.h \
    ../src/partyql/basicmath/Hadron.h \
    ../src/partyql/basicmath/Flat.h \
    ../src/partyql/basicmath/FinitePartition.h \
    ../src/partyql/basicmath/Descriptor.h \
    ../src/partyql/basicmath/Cone.h \
    ../src/partyql/basicmath/AffineTransformation.h \
    ../src/partyql/basicmath/factorization/HermiteNormalFormReduction.h \
    ../src/partyql/basicmath/factorization/EchelonReduction.h \
    ../src/partyql/basicmath/scalar/Rational.h \
    ../src/partyql/basicmath/scalar/Interval.h \
    ../src/partyql/basicmath/scalar/Integer.h \
    ../src/partyql/basicmath/scalar/ExtendedGCD.h \
    ../src/partyql/basicmath/scalar/Bool.h \
    ../src/partyql/cases/SCFG.h \
    ../src/partyql/cases/MatrixRecursion.h \
    ../src/partyql/core/Permutation.h \
    ../src/partyql/core/LinearProgramming.h \
    ../src/partyql/core/IntegerProgramming.h \
    ../src/partyql/core/HyperCube.h \
    ../src/partyql/core/HadronDescription.h \
    ../src/partyql/core/Hadron.h \
    ../src/partyql/core/Graph.h \
    ../src/partyql/core/Domain.h \
    ../src/partyql/core/Descriptor.h \
    ../src/partyql/core/Declarations.h \
    ../src/partyql/core/Config.h \
    ../src/partyql/core/Vektor.h \
    ../src/partyql/core/SparseVektor.h \
    ../src/partyql/core/SetRelation.h \
    ../src/partyql/core/scalar/Rational.h \
    ../src/partyql/core/scalar/Integer.h \
    ../src/partyql/core/scalar/ExtendedGCD.h \
    ../src/partyql/core/util/Util.h \
    ../src/partyql/core/util/Time.h \
    ../src/partyql/core/util/Printable.h \
    ../src/partyql/core/util/Html.h \
    ../src/partyql/core/util/DebugStream.h \
    ../src/partyql/graph/Graph.h \
    ../src/tests/Test_MultiListener.h \
    ../src/tests/Test_Funcs.h \
    ../src/utils/ToString.h \
    ../src/utils/SVGViewbox.h \
    ../src/utils/SVGEText.h \
    ../src/utils/SVGELine.h \
    ../src/utils/SVGElement.h \
    ../src/utils/SVGELattice.h \
    ../src/utils/SVGEGraph.h \
    ../src/utils/SVGEElementCollection.h \
    ../src/utils/SVGECircle.h \
    ../src/utils/SVG.h \
    ../src/utils/Point.h \
    ../src/utils/Html.h \
    ../src/utils/Color.h \
    ../cute/vstudio_listener.h \
    ../cute/ostream_listener.h \
    ../cute/ide_listener.h \
    ../cute/eclipse_listener.h \
    ../cute/cute.h \
    ../cute/cute_testmember.h \
    ../cute/cute_test.h \
    ../cute/cute_test_incarnate.h \
    ../cute/cute_suite.h \
    ../cute/cute_suite_test.h \
    ../cute/cute_listener.h \
    ../cute/cute_expect.h \
    ../cute/cute_equals.h \
    ../cute/cute_counting_listener.h \
    ../cute/cute_base.h \
    ../src/partyql/core/scalar/Bool.h

SOURCES += \
    ../src/Test.cpp \
    ../src/Render.cpp \
    ../src/partyql/Variable.cpp \
    ../src/partyql/StandardOperations.cpp \
    ../src/partyql/SpacePartitioner.cpp \
    ../src/partyql/Reference.cpp \
    ../src/partyql/RDMGraph.cpp \
    ../src/partyql/Navigator.cpp \
    ../src/partyql/CommunicationChannel.cpp \
    ../src/partyql/basicmath/Vektor.cpp \
    ../src/partyql/basicmath/VectorDescriptor.cpp \
    ../src/partyql/basicmath/Transitions.cpp \
    ../src/partyql/basicmath/SetRelation.cpp \
    ../src/partyql/basicmath/PolyhedralDomain.cpp \
    ../src/partyql/basicmath/Polyheder.cpp \
    ../src/partyql/basicmath/PartitionRelation.cpp \
    ../src/partyql/basicmath/Matrix.cpp \
    ../src/partyql/basicmath/IntegerProgramming.cpp \
    ../src/partyql/basicmath/HadronDescriptor.cpp \
    ../src/partyql/basicmath/Hadron.cpp \
    ../src/partyql/basicmath/Flat.cpp \
    ../src/partyql/basicmath/FinitePartition.cpp \
    ../src/partyql/basicmath/AffineTransformation.cpp \
    ../src/partyql/basicmath/factorization/HermiteNormalFormReduction.cpp \
    ../src/partyql/basicmath/factorization/EchelonReduction.cpp \
    ../src/partyql/basicmath/scalar/Rational.cpp \
    ../src/partyql/basicmath/scalar/Interval.cpp \
    ../src/partyql/basicmath/scalar/Integer.cpp \
    ../src/partyql/core/IntegerProgramming.cpp \
    ../src/partyql/core/HadronDescription.cpp \
    ../src/partyql/core/SetRelation.cpp \
    ../src/partyql/core/scalar/Rational.cpp \
    ../src/partyql/core/scalar/Integer.cpp \
    ../src/partyql/core/scalar/Bool.cpp \
    ../src/partyql/core/util/Util.cpp \
    ../src/partyql/core/util/Time.cpp \
    ../src/partyql/core/util/Html.cpp \
    ../src/partyql/core/util/DebugStream.cpp \
    ../src/partyql/graph/Vertex.cpp \
    ../src/partyql/graph/Graph.cpp \
    ../src/partyql/graph/Edge.cpp \
    ../src/tests/Test_Vector.cpp \
    ../src/tests/Test_StronglyConnectedComponents.cpp \
    ../src/tests/Test_SetRelation_NewHadron.cpp \
    ../src/tests/Test_RDMGraph.cpp \
    ../src/tests/Test_RADGraph.cpp \
    ../src/tests/Test_Polyheder.cpp \
    ../src/tests/Test_ModuleRelation.cpp \
    ../src/tests/Test_Module.cpp \
    ../src/tests/Test_Matrix.cpp \
    ../src/tests/Test_LatticeRelation.cpp \
    ../src/tests/Test_Lattice.cpp \
    ../src/tests/Test_Lattice_NewHadron.cpp \
    ../src/tests/Test_IntegerProgramming.cpp \
    ../src/tests/Test_Integer.cpp \
    ../src/tests/Test_Graph.cpp \
    ../src/tests/Test_Funcs.cpp \
    ../src/tests/Test_FourierMotzkin.cpp \
    ../src/tests/Test_FinitePartition.cpp \
    ../src/tests/Test_DeskriptorSet.cpp \
    ../src/tests/Test_ConnectedComponents.cpp \
    ../src/tests/Test_Cone.cpp \
    ../src/tests/Test_Computation.cpp \
    ../src/tests/Test_Computation_NewHadron.cpp \
    ../src/utils/ToString.cpp \
    ../src/utils/SVGElement.cpp \
    ../src/utils/SVGEGraph.cpp \
    ../src/utils/SVG.cpp \
    ../src/utils/Html.cpp \
    ../src/utils/Color.cpp

OTHER_FILES += \
    ../LICENSE-MIT.txt \
    ../cute/CUTE.txt































