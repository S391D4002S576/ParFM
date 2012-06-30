#include "../cute/cute.h"
#include "../cute/cute_suite.h"

#include "Test_Funcs.h"

//#include "partyql/graph/Graph.h"
//#include "partyql/basicmath/scalar/Integer.h"
//#include "partyql/basicmath/Vektor.h"
//#include "partyql/basicmath/Module.h"
//#include "partyql/basicmath/SetRelation.h"
//#include "partyql/basicmath/AffineTransformation.h"
//#include "partyql/basicmath/PolyhedralDomain.h"
//#include "../partyql/RDMGraph.h"
//#include "../partyql/HierarGraph.h"
#include "../partyql/Computation.h"
#include "../partyql/core/Hadron.h"
#include "../partyql/core/Descriptor.h"
#include "../partyql/core/Domain.h"
//#include "../partyql/Statement.h"

/*#include "../utils/SVG.h"
#include "../utils/SVGELattice.h"
#include "../utils/SVGViewbox.h"
#include "../utils/ToString.h"
#include "../utils/Html.h"*/

// #define ENABLE_PDF

using namespace AlgoTrans;

namespace Test_Computation_NewHadron {

	typedef CInteger I;
	typedef CVector<I> IVector;
	typedef CAffineTransformation<I> IAffineTransformation;
	typedef CModuleRelation<I>::T IModuleRelation;
	typedef CLatticeRelation<I>::T ILatticeRelation;

	typedef CComputation<I> IComputation;
	typedef CVariable<I> IVariable;
	typedef CStatement<I> IStatement;

	typedef CVertex<CPolyheder<I>, IAffineTransformation> RDMVertex;
	typedef CEdge<CPolyheder<I>, IAffineTransformation> RDMEdge;
	//typedef CRDMGraph<I> RDMGraph;

	typedef CVertex<CPolyheder<I>, ILatticeRelation> RADLatticeVertex;
	typedef CEdge<CPolyheder<I>, ILatticeRelation> RADLatticeEdge;
	typedef CRADGraph<I, ILatticeRelation> RADLatticeGraph;

	typedef CVertex<CPolyheder<I>, IModuleRelation> RADModuleVertex;
	typedef CEdge<CPolyheder<I>, IModuleRelation> RADModuleEdge;
	typedef CRADGraph<I, IModuleRelation> RADModuleGraph;

	//typedef CHierarNode<CComputation<I>, I> IHierarPart;
	typedef CHierarPart<CComputation<I>, I> IHierarPart;
	typedef CHierarLeafPart<CComputation<I>, I> IHierarLeafPart;

	RADModuleGraph* radgModule;
	RADLatticeGraph* radgLattice;

	vector<RADModuleVertex* > wModule;
	vector<RADLatticeVertex* > wLattice;

	vector<vector<vector<RADModuleEdge*> > > resultEdges;
	vector<vector<RADLatticeEdge*> > resultEdgesLattice;

	vector<ILattice*> distanceLattices;

	string verboseOutputPrefix, verboseOutputPath, fwInputPdf, initialLatticesPdf;
	CVector<int> lowerBounds, upperBounds;

	vector<IVariable*> v;

	IComputation setupComputation(int variables) {
		IComputation comp = IComputation();

		v.clear();
		for (int q = 0; q < variables; q++) v.push_back(&comp.addVariable());

		return comp;
	}

	/* 1D Sequential loop
	 * for ( i = 0 ;   ; i++) {
	 *     A[i] := A[i] + A[i - 1];
	 * }	 */
	void Trivial1DSequential_While() {
		IComputation comp = setupComputation(1);

		IStatement& s = comp.addStatement();

		comp.setOriginalExecutionOrder(IHierarPart::fromLexicoLinearOrderedSubNodes(comp, IPDC(false, 1, IV_(1)), 1,
				IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(2,  0,  1), IV_(2,  1,  0)), 1,
						new IHierarLeafPart(s)
				)
		));

		IAffineTransformation iat_i = IAT(1, IV_(2,  0,  1));
		IAffineTransformation iat_i_1 = IAT(1, IV_(2, -1,  1));
		comp.addReference(s, *v[0], WRITE, iat_i);
		comp.addReference(s, *v[0], READ,  iat_i);
		comp.addReference(s, *v[0], READ,  iat_i_1);

		// Assert we have the correct full iteration spaces
		ASSERT(s.getFullIterationSpace() == IPDC(false, 2, IV_(2,  0,  1), IV_(2,  1,  0)));

		// Perform partitioning
		DebugOutStream ofs =  DebugOutStream("Trivial1DSequential_While");
		comp.osHtml = &ofs;

		comp.computeImplementation();

		IHierarScatterPart* timeScatterPart = dynamic_cast<IHierarScatterPart* >(comp.executionOrder);
		ASSERT(timeScatterPart != NULL);
		ASSERT(timeScatterPart->isTimePart());
		ASSERT_EQ(timeScatterPart->getScatterData(), IH(G, true, 1, IV_(2, 1, 0)) + IH(G, false, 1, IV_(2, 0, 1)));

		IHierarLeafPart* leafPart = dynamic_cast<IHierarLeafPart* >(timeScatterPart->getChildPart());
		ASSERT(leafPart != NULL);
		ASSERT(&leafPart->getStatement() == &s);	}

	/* 1D Sequential loop
	 * for ( i = 0 ;   ; i++) {
	 *     A[i] := A[i] + A[i - 2];
	 * }	 */
	void Trivial1DSequential_While_L() {
		IComputation comp = setupComputation(1);

		IStatement& s = comp.addStatement();

		comp.setOriginalExecutionOrder(IHierarPart::fromLexicoLinearOrderedSubNodes(comp, IPDC(false, 0, 0), 1,
				IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 1, IV_(2,  0,  1)), 1,
						new IHierarLeafPart(s)
				)
		));

		IAffineTransformation iat_i = IAT(1, IV_(2,  0,  1));
		IAffineTransformation iat_i_2 = IAT(1, IV_(2, -2,  1));
		comp.addReference(s, *v[0], WRITE, iat_i);
		comp.addReference(s, *v[0], READ,  iat_i);
		comp.addReference(s, *v[0], READ,  iat_i_2);

		// Assert we have the correct full iteration spaces
		ASSERT(s.getFullIterationSpace() == IPDC(false, 1, IV_(2,  0,  1)));

		// Perform partitioning
		DebugOutStream ofs =  DebugOutStream("Trivial1DSequential_While_L");
		comp.osHtml = &ofs;

		comp.computeImplementation();

		ASSERT(*comp.spaceSolutions == IH(G, true, 1, IV_(2, 1, 0)));
		ASSERT(*comp.timeSolutions == IH(G, true, 1, IV_(2, 1, 0)) + IH(G, false, 1, IV_(2, 0, 1)));
	}

	/* 1D Sequential loop
	 * with (n >= 1):
	 * for ( i = 0 ; i <= n ; i++) {
	 *     A[i] := A[i] + A[i - 1];
	 * }
	 */
	void Trivial1DSequential() {
		IComputation comp = setupComputation(1);

		IStatement& s = comp.addStatement();
		comp.setOriginalExecutionOrder(IHierarPart::fromLexicoLinearOrderedSubNodes(comp, IPDC(false, 2, IV_(2, -1,  1), IV_(2, 1,  0)), 1,
				IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(3,  0,  0,  1),
                        										IV_(3,  0,  1, -1)), 1,
						new IHierarLeafPart(s)
				)
		));

		IAffineTransformation iat_i = IAT(1, IV_(3,  0,  0,  1));
		IAffineTransformation iat_i_1 = IAT(1, IV_(3, -1,  0,  1));
		comp.addReference(s, *v[0], WRITE, iat_i);
		comp.addReference(s, *v[0], READ,  iat_i);
		comp.addReference(s, *v[0], READ,  iat_i_1);

		// Assert we have the correct full iteration spaces
		ASSERT_EQ(s.getFullIterationSpace(), IPDC(false, 4, IV_(3, -1,  1,  0),
 	                                                        IV_(3,  1,  0,  0),
											 	            IV_(3,  0,  0,  1),
												            IV_(3,  0,  1, -1)));

		// Perform partitioning
		DebugOutStream ofs =  DebugOutStream("Trivial1DSequential");
		comp.osHtml = &ofs;

		comp.computeImplementation();
		//comp.performAffineSetSpacePartitioning();
		//comp.performAffineLatticeSpacePartitioning();

		ASSERT_EQ(*comp.spaceSolutions, IH(G, true, 2, IV_(3, 1, 0, 0), IV_(3, 0, 1, 0)));
		ASSERT_EQ(*comp.timeSolutions, IH(G, true, 2, IV_(3, 1, 0, 0), IV_(3, 0, 1, 0)) + IH(G, false, 1, IV_(3, 0, 0, 1)));
	}

	/* for i := 0 to 77 do begin
	 *   for j := 0 to 77 do begin
	 *     A[i][j] := A[i - 1][j] + A[i][j - 1]
	 *   end;
	 * end;
	 */
	void simple2DUniform() {
		IComputation comp = setupComputation(1);

		IStatement& s = comp.addStatement();

		comp.setOriginalExecutionOrder(IHierarPart::fromLexicoLinearOrderedSubNodes(comp, IPDC(false, 1, IV_(1)), 1,
			IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(2,  0,  1),
                      										IV_(2, 77, -1)), 1,
				IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(3,  0,  0,  1),
                        										IV_(3, 77,  0, -1)), 1,
						new IHierarLeafPart(s)
				)
			)
		));

		comp.addReference(s, *v[0], WRITE, IAT(2, IV_(3,  0,  1,  0),
				                                  IV_(3,  0,  0,  1)));
		comp.addReference(s, *v[0], READ,  IAT(2, IV_(3, -1,  1,  0),
				                                  IV_(3,  0,  0,  1)));
		comp.addReference(s, *v[0], READ,  IAT(2, IV_(3,  0,  1,  0),
				                                  IV_(3, -1,  0,  1)));

		// Perform partitioning
		DebugOutStream ofs =  DebugOutStream("simple2DUniform");
		comp.osHtml = &ofs;

		comp.computeImplementation();

		ASSERT_EQ(*comp.spaceSolutions, IH(G, true, 1, IV_(3, 1, 0, 0)));
		ASSERT_EQ(*comp.timeSolutions, IH(G, true, 1, IV_(3, 1, 0, 0)) + IH(G, false, 2, IV_(3, 0, 1, 0), IV_(3, 0, 0, 1)));
	}

	/* for i >= 0  do begin
	 *   for j >= 0 do begin
	 *     A[i][j] := A[i - 1][j] + A[i][j - 1]
	 *   end;
	 * end;
	 */
	void simple2DUniform_Unbounded() {
		IComputation comp = setupComputation(1);

		IStatement& s = comp.addStatement();

		comp.setOriginalExecutionOrder(IHierarPart::fromLexicoLinearOrderedSubNodes(comp, IPDC(false, 1, IV_(1)), 1,
			IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 1, IV_(2,  0,  1)), 1,
				IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 1, IV_(3,  0,  0,  1)), 1,
						new IHierarLeafPart(s)
				)
			)
		));

		comp.addReference(s, *v[0], WRITE, IAT(2, IV_(3,  0,  1,  0),
				                                  IV_(3,  0,  0,  1)));
		comp.addReference(s, *v[0], READ,  IAT(2, IV_(3, -1,  1,  0),
				                                  IV_(3,  0,  0,  1)));
		comp.addReference(s, *v[0], READ,  IAT(2, IV_(3,  0,  1,  0),
				                                  IV_(3, -1,  0,  1)));

		// Perform partitioning
		DebugOutStream ofs =  DebugOutStream("simple2DUniform_Unbounded");
		comp.osHtml = &ofs;

		comp.computeImplementation();

		ASSERT_EQ(*comp.spaceSolutions, IH(G, true, 1, IV_(3, 1, 0, 0)));
		ASSERT_EQ(*comp.timeSolutions, IH(G, true, 1, IV_(3, 1, 0, 0)) + IH(G, false, 2, IV_(3, 0, 1, 0), IV_(3, 0, 0, 1)));
	}

	/* Trivial 1D parallelism
	 * for ( i = 0 ; i <= n ; i++) {
	 *     A[i] := A[i];
	 * }
	 */
	void LexiCone() {
		IComputation comp = setupComputation(1);

		CDisjunctiveDomain<Polyheder>  pd = comp.getLexiCone(0, 1, 1, 1, true);

		CDisjunctiveDomain<Polyheder>  exp = IPDC(false, 2, IV_(3, 0, -1, 1), IV_(3, 1, 0, 0));

		ASSERT_EQ(pd, exp);
	}

	/* Trivial 1D parallelism
	 * for ( i = 0 ; i <= n ; i++) {
	 *     A[i] := A[i];
	 * }
	 */
	void Trivial1DParallelism() {
		IComputation comp = setupComputation(1);

		IStatement& sC = comp.addStatement();

		comp.setOriginalExecutionOrder(IHierarPart::fromLexicoLinearOrderedSubNodes(comp,
			IPDC(false, 2, IV_(2, 0,  1),   // n >= 0
                           IV_(2, 1,  0)), 1,
			IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(3,  0,  0,  1),   // i >= 0
                    										            IV_(3,  0,  1, -1)), 1,
						new IHierarLeafPart(sC)
			)
		));

		IAffineTransformation iat_i_j = IAT(1, IV_(3,  0,  0,  1));
		comp.addReference(sC, *v[0], WRITE, iat_i_j);
		comp.addReference(sC, *v[0], READ,  iat_i_j);

		// Assert we have the correct full iteration spaces
		ASSERT_EQ(sC.getFullIterationSpace(), IPDC(false, 3, IV_(3,  1,  0,  0),
												             IV_(3,  0,  0,  1),
												             IV_(3,  0,  1, -1)));

		// Perform partitioning
		DebugOutStream ofs =  DebugOutStream("Trivial1DParallelism");
		comp.osHtml = &ofs;

		comp.computeImplementation();

		ASSERT_EQ(*comp.spaceSolutions, IH(G, true, 3, IV_(3, 1, 0, 0),
				                                       IV_(3, 0, 1, 0),
													   IV_(3, 0, 0, 1)));

		ASSERT_EQ(*comp.timeSolutions, IH(G, true, 3, IV_(3, 1, 0, 0),
				                                      IV_(3, 0, 1, 0),
												      IV_(3, 0, 0, 1)));
	}

	/* Trivial 1D parallelism
	 * for ( i = 0 ;  ; i++) {
	 *     A[i] := A[i] + B[i];
	 *     B[i] := B[i] + A[i];
	 * }
	 */
	void Trivial1DParallelism_2S() {
		IComputation comp = setupComputation(2);

		IStatement& sC = comp.addStatement();
		IStatement& sD = comp.addStatement();

		comp.setOriginalExecutionOrder(IHierarPart::fromLexicoLinearOrderedSubNodes(comp, IPDC(false, 1, IV_(1,  1)), 1,
			IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 1, IV_(2,  0,  1)), 2,
						new IHierarLeafPart(sC),
						new IHierarLeafPart(sD)
			)
		));

		IAffineTransformation iat_i_j = IAT(1, IV_(2,  0,  1));
		comp.addReference(sC, *v[0], WRITE, iat_i_j);
		comp.addReference(sC, *v[0], READ,  iat_i_j);
		comp.addReference(sC, *v[1], READ, iat_i_j);
		comp.addReference(sD, *v[1], WRITE, iat_i_j);
		comp.addReference(sD, *v[1], READ,  iat_i_j);
		comp.addReference(sD, *v[0], READ,  iat_i_j);

		// Assert we have the correct full iteration spaces
		ASSERT_EQ(sC.getFullIterationSpace(), IPDC(false, 2, IV_(2,  1,  0),
												             IV_(2,  0,  1)));
		ASSERT_EQ(sD.getFullIterationSpace(), IPDC(false, 2, IV_(2,  1,  0),
												             IV_(2,  0,  1)));

		// Perform partitioning
		DebugOutStream ofs =  DebugOutStream("Trivial1DParallelism_2S");
		comp.osHtml = &ofs;

		comp.computeImplementation();

		ASSERT_EQ(*comp.spaceSolutions, IH(G, true, 2, IV_(4, 1, 0, 1, 0),
				                                       IV_(4, 0, 1, 0, 1)));

		ASSERT_EQ(*comp.timeSolutions, IH(G, true, 2, IV_(4, 1, 0, 1, 0),
                                                      IV_(4, 0, 1, 0, 1))
                                       + IH(G, false, 2, IV_(4, 0, 0, 1, 0),
                                                      IV_(4, 0, 0, 0, 1)));
	}

	/* Trivial 1D parallelism
	 * for ( i = 0 ; i <= n ; i++) {
	 *     A[i] := A[i] + B[i];
	 *     B[i] := B[i] + A[i];
	 * }
	 */
	void Trivial1DParallelism_2S_1P() {
		IComputation comp = setupComputation(2);

		IStatement& sC = comp.addStatement();
		IStatement& sD = comp.addStatement();

		comp.setOriginalExecutionOrder(IHierarPart::fromLexicoLinearOrderedSubNodes(comp,
			IPDC(false, 2, IV_(2, 0,  1),   // n >= 0
                           IV_(2, 1,  0)), 1,
			IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(3,  0,  0,  1),   // i >= 0
                    										            IV_(3,  0,  1, -1)), 2,
						new IHierarLeafPart(sC),
						new IHierarLeafPart(sD)
			)
		));

		IAffineTransformation iat_i_j = IAT(1, IV_(3,  0,  0,  1));
		comp.addReference(sC, *v[0], WRITE, iat_i_j);
		comp.addReference(sC, *v[0], READ,  iat_i_j);
		comp.addReference(sC, *v[1], READ, iat_i_j);
		comp.addReference(sD, *v[1], WRITE, iat_i_j);
		comp.addReference(sD, *v[1], READ,  iat_i_j);
		comp.addReference(sD, *v[0], READ,  iat_i_j);

		// Assert we have the correct full iteration spaces
		ASSERT_EQ(sC.getFullIterationSpace(), IPDC(false, 3, IV_(3,  1,  0,  0),
												             IV_(3,  0,  0,  1),
												             IV_(3,  0,  1, -1)));
		ASSERT_EQ(sD.getFullIterationSpace(), IPDC(false, 3, IV_(3,  1,  0,  0),
												             IV_(3,  0,  0,  1),
												             IV_(3,  0,  1, -1)));

		// Perform partitioning
		DebugOutStream ofs =  DebugOutStream("Trivial1DParallelism_2S_1P");
		comp.osHtml = &ofs;

		comp.computeImplementation();

		ASSERT_EQ(*comp.spaceSolutions, IH(G, true, 3, IV_(6, 1, 0, 0, 1, 0, 0),
				                                       IV_(6, 0, 1, 0, 0, 1, 0),
													   IV_(6, 0, 0, 1, 0, 0, 1)));

		ASSERT_EQ(*comp.timeSolutions, IH(G, true, 3, IV_(6, 1, 0, 0, 1, 0, 0),
				                                      IV_(6, 0, 1, 0, 0, 1, 0),
												      IV_(6, 0, 0, 1, 0, 0, 1))
				                       + IH(G, false, 3, IV_(6, 0, 0, 0, 1, 0, 0),
												        IV_(6, 0, 0, 1, 0, 1, 0),
												        IV_(6, 0, 0, 0, 0, 0, 1)
												      ));
	}

	/* Trivial 1D parallelism
	 * for ( i = 0 ; i <= n ; i++) {
	 *     A[i] := A[i] + B[i];
	 *     B[i] := B[i] + A[i];
	 * }
	 */
	void Trivial1DParallelism_2S_1P_B() {
		IComputation comp = setupComputation(2);

		IStatement& sC = comp.addStatement();
		IStatement& sD = comp.addStatement();

		comp.setOriginalExecutionOrder(IHierarPart::fromLexicoLinearOrderedSubNodes(comp,
			IPDC(false, 1, IV_(2, 1,  0)), 1,
			IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(3,  0,  0,  1),   // i >= 0
                    										            IV_(3,  0,  1, -1)), 2,
						new IHierarLeafPart(sC),
						new IHierarLeafPart(sD)
			)
		));

		IAffineTransformation iat_i_j = IAT(1, IV_(3,  0,  0,  1));
		comp.addReference(sC, *v[0], WRITE, iat_i_j);
		comp.addReference(sC, *v[0], READ,  iat_i_j);
		comp.addReference(sC, *v[1], READ, iat_i_j);
		comp.addReference(sD, *v[1], WRITE, iat_i_j);
		comp.addReference(sD, *v[1], READ,  iat_i_j);
		comp.addReference(sD, *v[0], READ,  iat_i_j);

		// Assert we have the correct full iteration spaces
		ASSERT_EQ(sC.getFullIterationSpace(), IPDC(false, 3, IV_(3,  1,  0,  0),
												             IV_(3,  0,  0,  1),
												             IV_(3,  0,  1, -1)));
		ASSERT_EQ(sD.getFullIterationSpace(), IPDC(false, 3, IV_(3,  1,  0,  0),
												             IV_(3,  0,  0,  1),
												             IV_(3,  0,  1, -1)));

		// Perform partitioning
		DebugOutStream ofs =  DebugOutStream("Trivial1DParallelism_2S_1P_B");
		comp.osHtml = &ofs;

		comp.computeImplementation();

		IHierarScatterPart* spaceScatterPart = dynamic_cast<IHierarScatterPart* >(comp.executionOrder);
		ASSERT(!spaceScatterPart->isTimePart());
		ASSERT(spaceScatterPart != NULL);
		ASSERT_EQ(spaceScatterPart->getScatterData(), IH(G, true, 3, IV_(6, 1, 0, 0, 1, 0, 0),
				                                       IV_(6, 0, 1, 0, 0, 1, 0),
													   IV_(6, 0, 0, 1, 0, 0, 1)));

		IHierarFinitePart* spaceFinitePart = dynamic_cast<IHierarFinitePart* >(spaceScatterPart->getChildPart());
		ASSERT(spaceFinitePart != NULL);
		ASSERT_EQ(I(spaceFinitePart->getChildCount()), I(2));

		vector<CVertex<IHierarPart*, Bool>* > childParts = spaceFinitePart->cellGraph.getVertices();

		ASSERT(dynamic_cast<IHierarLeafPart* >(childParts[0]->getData()) != NULL);
		ASSERT(dynamic_cast<IHierarLeafPart* >(childParts[1]->getData()) != NULL);
/*		CHierarScatterPart<CComputation<CInteger>, CInteger >* timeScatterPart
		= dynamic_cast<CHierarScatterPart<CComputation<CInteger>, CInteger>* >(comp.executionOrder);*/

	}

	/* Trivial 1D parallelism
	 * for ( i = 0 ; i <= n ; i++) {
	 *     A[i] := A[i] + A[i];
	 *     B[i] := B[i] + B[i];
	 * }
	 */
	void Trivial1DParallelism_2S_1P_C() {
		IComputation comp = setupComputation(2);

		IStatement& sC = comp.addStatement();
		IStatement& sD = comp.addStatement();

		comp.setOriginalExecutionOrder(IHierarPart::fromLexicoLinearOrderedSubNodes(comp,
			IPDC(false, 1, IV_(2, 1,  0)), 1,
			IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(3,  0,  0,  1),   // i >= 0
                    										            IV_(3,  0,  1, -1)), 2,
						new IHierarLeafPart(sC),
						new IHierarLeafPart(sD)
			)
		));

		IAffineTransformation iat_i_j = IAT(1, IV_(3,  0,  0,  1));
		comp.addReference(sC, *v[0], WRITE, iat_i_j);
		comp.addReference(sC, *v[0], READ,  iat_i_j);
		comp.addReference(sC, *v[0], READ, iat_i_j);
		comp.addReference(sD, *v[1], WRITE, iat_i_j);
		comp.addReference(sD, *v[1], READ,  iat_i_j);
		comp.addReference(sD, *v[1], READ,  iat_i_j);

		// Assert we have the correct full iteration spaces
		ASSERT_EQ(sC.getFullIterationSpace(), IPDC(false, 3, IV_(3,  1,  0,  0),
												             IV_(3,  0,  0,  1),
												             IV_(3,  0,  1, -1)));
		ASSERT_EQ(sD.getFullIterationSpace(), IPDC(false, 3, IV_(3,  1,  0,  0),
												             IV_(3,  0,  0,  1),
												             IV_(3,  0,  1, -1)));

		// Perform partitioning
		DebugOutStream ofs =  DebugOutStream("Trivial1DParallelism_2S_1P");
		comp.osHtml = &ofs;

		comp.computeImplementation();

		IHierarFinitePart* finPart = dynamic_cast<IHierarFinitePart*>(comp.executionOrder);
		ASSERT_EQ(I(finPart->getChildCount()), I(2));

		vector<CVertex<IHierarPart*, Bool>* > childParts = finPart->cellGraph.getVertices();

		ASSERT(dynamic_cast<IHierarScatterPart* >(childParts[0]->getData()) != NULL);
		IHierarScatterPart* cA = dynamic_cast<IHierarScatterPart* >(childParts[0]->getData());
		ASSERT(!cA->isTimePart());
		ASSERT(dynamic_cast<IHierarScatterPart* >(childParts[1]->getData()) != NULL);
		IHierarScatterPart* cB = dynamic_cast<IHierarScatterPart* >(childParts[1]->getData());
		ASSERT(!cB->isTimePart());

		ASSERT_EQ(cA->getScatterData(), IH(G, true, 3, IV_(3, 1, 0, 0),
				                                       IV_(3, 0, 1, 0),
													   IV_(3, 0, 0, 1)));
		ASSERT_EQ(cB->getScatterData(), IH(G, true, 3, IV_(3, 1, 0, 0),
				                                       IV_(3, 0, 1, 0),
													   IV_(3, 0, 0, 1)));

		ASSERT(dynamic_cast<IHierarLeafPart* >(cA->getChildPart()) != NULL);
		IHierarLeafPart* lA = dynamic_cast<IHierarLeafPart* >(cA->getChildPart());
		ASSERT(dynamic_cast<IHierarLeafPart* >(cB->getChildPart()) != NULL);
		IHierarLeafPart* lB = dynamic_cast<IHierarLeafPart* >(cB->getChildPart());

		ASSERT((&lA->getStatement() == &sC) || (&lA->getStatement() == &sD));
		ASSERT(((&lB->getStatement() == &sC) || (&lB->getStatement() == &sD)) && (&lA->getStatement() != &lB->getStatement())) ;
	}

	/* Trivial 2D parallelism
	 * for ( i = 0 ; i <= n ; i++) {
	 * 	 for ( j = 0 ; j <= m ; j++) {
	 *     A[i][j] := A[i][j] + B[i][j];
	 *   }
	 * }
	 */
	void Trivial2DParallelism_SingleStatement() {
		IComputation comp = setupComputation(2);

		IStatement& sC = comp.addStatement();

		comp.setOriginalExecutionOrder(IHierarPart::fromLexicoLinearOrderedSubNodes(comp,
			IPDC(false, 3, IV_(3,  0,  1,  0),   // n >= 0
                           IV_(3,  0,  0,  1),
                           IV_(3,  1,  0,  0)), 1,
			IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(4,  0,  0,  0,  1),   // i >= 0
                    										            IV_(4,  0,  1,  0, -1)), 1,
				IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(5,  0,  0,  0,  0,  1),   // j >= 0
                        										            IV_(5,  0,  0,  1,  0, -1)), 1,
						new IHierarLeafPart(sC)
				)
			)
		));

		IAffineTransformation iat_i_j = IAT(2, IV_(5,  0,  0,  0,  1,  0),
                							   IV_(5,  0,  0,  0,  0,  1));
		comp.addReference(sC, *v[0], WRITE, iat_i_j);
		comp.addReference(sC, *v[0], READ,  iat_i_j);
		comp.addReference(sC, *v[0], READ,  iat_i_j);

		// Assert we have the correct full iteration spaces
		ASSERT_EQ(sC.getFullIterationSpace(), IPDC(false, 7, IV_(5, 0,  1,  0,  0,  0),
												           IV_(5, 0,  0,  1,  0,  0),
												           IV_(5,  0,  0,  0,  1,  0),
												           IV_(5,  0,  1,  0, -1,  0),
												           IV_(5,  0,  0,  0,  0,  1),
												           IV_(5,  1,  0,  0,  0,  0),
												           IV_(5,  0,  0,  1,  0, -1)));

		// Perform partitioning
		DebugOutStream ofs =  DebugOutStream("Trivial2DParallelism_SingleStatement");
		comp.osHtml = &ofs;

		comp.computeImplementation();

		ASSERT_EQ(*comp.spaceSolutions, IH(G, true, 5, IV_(5, 1, 0, 0, 0, 0),
				                                       IV_(5, 0, 1, 0, 0, 0),
													   IV_(5, 0, 0, 1, 0, 0),
				                                       IV_(5, 0, 0, 0, 1, 0),
				                                       IV_(5, 0, 0, 0, 0, 1)));

		ASSERT_EQ(*comp.timeSolutions, IH(G, true, 5, IV_(5, 1, 0, 0, 0, 0),
				                                      IV_(5, 0, 1, 0, 0, 0),
													  IV_(5, 0, 0, 1, 0, 0),
				                                      IV_(5, 0, 0, 0, 1, 0),
				                                      IV_(5, 0, 0, 0, 0, 1)));
		//comp.performAffineLatticeSpacePartitioning();
	}

	/* Trivial 2D parallelism
	 * for ( i = 0 ; i <= n ; i++) {
	 * 	 for ( j = 0 ; j <= m ; j++) {
	 *     A[i][j] := A[i][j] + B[i][j];
	 *     B[i][j] := A[i][j] * B[i][j];
	 *   }
	 * }
	 */
	void Trivial2DParallelism() {
		IComputation comp = setupComputation(2);

		IStatement& sC = comp.addStatement();
		IStatement& sD = comp.addStatement();

		comp.setOriginalExecutionOrder(IHierarPart::fromLexicoLinearOrderedSubNodes(comp,
			IPDC(false, 3, IV_(3, 0,  1,  0),   // n >= 1
                     IV_(3, 0,  0,  1),
                     IV_(3,  1,  0,  0)), 1,
			IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(4,  0,  0,  0,  1),   // i >= 0
                    										IV_(4,  0,  1,  0, -1)), 1,
				IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(5,  0,  0,  0,  0,  1),   // j >= 0
                        										IV_(5,  0,  0,  1,  0, -1)), 2,
						new IHierarLeafPart(sC),
						new IHierarLeafPart(sD)
				)
			)
		));

		IAffineTransformation iat_i_j = IAT(2, IV_(5,  0,  0,  0,  1,  0),
                							   IV_(5,  0,  0,  0,  0,  1));
		comp.addReference(sC, *v[0], WRITE, iat_i_j);
		comp.addReference(sC, *v[0], READ,  iat_i_j);
		comp.addReference(sC, *v[0], READ,  iat_i_j);
		comp.addReference(sD, *v[0], WRITE, iat_i_j);
		comp.addReference(sD, *v[0], READ,  iat_i_j);
		comp.addReference(sD, *v[0], READ,  iat_i_j);

		// Assert we have the correct full iteration spaces
		ASSERT_EQ(sC.getFullIterationSpace(), IPDC(false, 7, IV_(5, 0,  1,  0,  0,  0),
												           IV_(5, 0,  0,  1,  0,  0),
												           IV_(5,  0,  0,  0,  1,  0),
												           IV_(5,  0,  1,  0, -1,  0),
												           IV_(5,  0,  0,  0,  0,  1),
												           IV_(5,  1,  0,  0,  0,  0),
												           IV_(5,  0,  0,  1,  0, -1)));

		ASSERT_EQ(sD.getFullIterationSpace(), IPDC(false, 7, IV_(5,  0,  1,  0,  0,  0),
														   IV_(5,  0,  0,  1,  0,  0),
														   IV_(5,  0,  0,  0,  1,  0),
														   IV_(5,  0,  1,  0, -1,  0),
														   IV_(5,  0,  0,  0,  0,  1),
												           IV_(5,  1,  0,  0,  0,  0),
														   IV_(5,  0,  0,  1,  0, -1)));

		// Perform partitioning
		DebugOutStream ofs =  DebugOutStream("Trivial2DParallelism");
		comp.osHtml = &ofs;

		comp.computeImplementation();

		ASSERT_EQ(*comp.spaceSolutions, IH(G, true, 5, IV_(10, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0),
				                                      IV_(10, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0),
													  IV_(10, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0),
				                                      IV_(10, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0),
				                                      IV_(10, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1)));

		ASSERT_EQ(*comp.timeSolutions, IH(G, true, 5, IV_(10, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0),
				                                      IV_(10, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0),
													  IV_(10, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0),
				                                      IV_(10, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0),
				                                      IV_(10, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1))
				                       + IH(G, false, 5,
				                                      IV_(10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1),
				                                      IV_(10, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0),
				                                      IV_(10, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0),
				                                      IV_(10, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0),
				                                      IV_(10, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0)
				                                      ));
	}

	/* Trivial 2D parallelism
	 * for ( i = 0 ; i <= n ; i++) {
	 * 	 for ( j = 0 ; j <= m ; j++) {
	 *     A[i][j] := A[i][j] + B[i][j];
	 *     B[i][j] := A[i][j] * B[i][j];
	 *   }
	 * }
	 */
	void Trivial2DParallelism_B() {
		IComputation comp = setupComputation(2);

		IStatement& sC = comp.addStatement();
		IStatement& sD = comp.addStatement();

		comp.setOriginalExecutionOrder(IHierarPart::fromLexicoLinearOrderedSubNodes(comp,
			IPDC(false, 1, IV_(3,  1,  0,  0)), 1,
			IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(4,  0,  0,  0,  1),   // i >= 0
                    										IV_(4,  0,  1,  0, -1)), 1,
				IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(5,  0,  0,  0,  0,  1),   // j >= 0
                        										IV_(5,  0,  0,  1,  0, -1)), 2,
						new IHierarLeafPart(sC),
						new IHierarLeafPart(sD)
				)
			)
		));

		IAffineTransformation iat_i_j = IAT(2, IV_(5,  0,  0,  0,  1,  0),
                							   IV_(5,  0,  0,  0,  0,  1));
		comp.addReference(sC, *v[0], WRITE, iat_i_j);
		comp.addReference(sC, *v[0], READ,  iat_i_j);
		comp.addReference(sC, *v[0], READ,  iat_i_j);
		comp.addReference(sD, *v[0], WRITE, iat_i_j);
		comp.addReference(sD, *v[0], READ,  iat_i_j);
		comp.addReference(sD, *v[0], READ,  iat_i_j);

		// Assert we have the correct full iteration spaces
		ASSERT_EQ(sC.getFullIterationSpace(), IPDC(false, 7, IV_(5, 0,  1,  0,  0,  0),
												           IV_(5, 0,  0,  1,  0,  0),
												           IV_(5,  0,  0,  0,  1,  0),
												           IV_(5,  0,  1,  0, -1,  0),
												           IV_(5,  0,  0,  0,  0,  1),
												           IV_(5,  1,  0,  0,  0,  0),
												           IV_(5,  0,  0,  1,  0, -1)));

		ASSERT_EQ(sD.getFullIterationSpace(), IPDC(false, 7, IV_(5,  0,  1,  0,  0,  0),
														   IV_(5,  0,  0,  1,  0,  0),
														   IV_(5,  0,  0,  0,  1,  0),
														   IV_(5,  0,  1,  0, -1,  0),
														   IV_(5,  0,  0,  0,  0,  1),
												           IV_(5,  1,  0,  0,  0,  0),
														   IV_(5,  0,  0,  1,  0, -1)));

		// Perform partitioning
		DebugOutStream ofs =  DebugOutStream("Trivial2DParallelism");
		comp.osHtml = &ofs;

		comp.computeImplementation();

		ASSERT_EQ(*comp.spaceSolutions, IH(G, true, 5, IV_(10, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0),
				                                      IV_(10, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0),
													  IV_(10, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0),
				                                      IV_(10, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0),
				                                      IV_(10, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1)));

		ASSERT_EQ(*comp.timeSolutions, IH(G, true, 5, IV_(10, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0),
				                                      IV_(10, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0),
													  IV_(10, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0),
				                                      IV_(10, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0),
				                                      IV_(10, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1))
				                       + IH(G, false, 5,
				                                      IV_(10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1),
				                                      IV_(10, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0),
				                                      IV_(10, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0),
				                                      IV_(10, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0),
				                                      IV_(10, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0)
				                                      ));
		//comp.performAffineLatticeSpacePartitioning();
	}

	/* Trivial 2D parallelism
	 * with (n >= 1, m >= 2):
	 * for ( i = 0 ; i <= n ; i++) {
	 * 	 for ( j = 0 ; j <= m ; j++) {
	 *     A[i][j] := A[i][j] + B[i][j];
	 *     B[i][j] := A[i][j] * B[i][j];
	 *   }
	 * }
	 */
	void Trivial2DParallelism_C() {
		IComputation comp = setupComputation(2);

		IStatement& sC = comp.addStatement();
		IStatement& sD = comp.addStatement();

		comp.setOriginalExecutionOrder(IHierarPart::fromLexicoLinearOrderedSubNodes(comp,
			IPDC(false, 3, IV_(3, -1,  1,  0),   // n >= 1
                     IV_(3, -2,  0,  1),
                     IV_(3,  1,  0,  0)), 1,
			IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(4,  0,  0,  0,  1),   // i >= 0
                    										IV_(4,  0,  1,  0, -1)), 1,
				IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(5,  0,  0,  0,  0,  1),   // j >= 0
                        										IV_(5,  0,  0,  1,  0, -1)), 2,
						new IHierarLeafPart(sC),
						new IHierarLeafPart(sD)
				)
			)
		));

		IAffineTransformation iat_i_j = IAT(2, IV_(5,  0,  0,  0,  1,  0),
                							   IV_(5,  0,  0,  0,  0,  1));
		comp.addReference(sC, *v[0], WRITE, iat_i_j);
		comp.addReference(sC, *v[0], READ,  iat_i_j);
		comp.addReference(sC, *v[0], READ,  iat_i_j);
		comp.addReference(sD, *v[0], WRITE, iat_i_j);
		comp.addReference(sD, *v[0], READ,  iat_i_j);
		comp.addReference(sD, *v[0], READ,  iat_i_j);

		// Assert we have the correct full iteration spaces
		ASSERT_EQ(sC.getFullIterationSpace(), IPDC(false, 7, IV_(5, -1,  1,  0,  0,  0),
												           IV_(5, -2,  0,  1,  0,  0),
												           IV_(5,  0,  0,  0,  1,  0),
												           IV_(5,  0,  1,  0, -1,  0),
												           IV_(5,  0,  0,  0,  0,  1),
												           IV_(5,  1,  0,  0,  0,  0),
												           IV_(5,  0,  0,  1,  0, -1)));

		ASSERT_EQ(sD.getFullIterationSpace(), IPDC(false, 7, IV_(5, -1,  1,  0,  0,  0),
														   IV_(5, -2,  0,  1,  0,  0),
														   IV_(5,  0,  0,  0,  1,  0),
														   IV_(5,  0,  1,  0, -1,  0),
														   IV_(5,  0,  0,  0,  0,  1),
												           IV_(5,  1,  0,  0,  0,  0),
														   IV_(5,  0,  0,  1,  0, -1)));

		// Perform partitioning
		DebugOutStream ofs =  DebugOutStream("Trivial2DParallelism");
		comp.osHtml = &ofs;

		comp.computeImplementation();

		ASSERT(*comp.spaceSolutions == IH(G, true, 5, IV_(10, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0),
				                                      IV_(10, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0),
													  IV_(10, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0),
				                                      IV_(10, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0),
				                                      IV_(10, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1)));

		ASSERT_EQ(*comp.timeSolutions, IH(G, true, 5, IV_(10, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0),
				                                      IV_(10, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0),
													  IV_(10, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0),
				                                      IV_(10, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0),
				                                      IV_(10, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1))
				                       + IH(G, false, 5,
				                                      IV_(10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1),
				                                      IV_(10, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0),
				                                      IV_(10, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0),
				                                      IV_(10, 0, 0, 0, 0, 0, -2, 1, 0, 0, 0),
				                                      IV_(10, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0),
				                                      IV_(10, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0),
				                                      IV_(10, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0)
				                                      ));
	}

	/* From: Partitioning Loops with Variable Dependence Distances (Yu, D'Hollander)
	 * for i := -9 to 9 do begin
	 *   for j := -9 to 9 do begin
	 *     A[3i + 1][2i + j - 1] := A[i + 3][j + 1]^2
	 *   end;
	 * end;
	 */
	void Yu_1() {
		IComputation comp = setupComputation(1);

		IStatement& s = comp.addStatement();

		comp.addReference(s, *v[0], WRITE, IAT(2, IV_(3,  1,  3,  0),
				                                  IV_(3, -1,  2,  1)));
		comp.addReference(s, *v[0], READ,  IAT(2, IV_(3,  3,  1,  0),
				                                  IV_(3,  1,  0,  1)));

		comp.setOriginalExecutionOrder(IHierarPart::fromLexicoLinearOrderedSubNodes(comp, IPDC(false, 0, 0), 1,
			IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(2,  9,  1),   // i >= -9
                    										IV_(2, -9, -1)), 1,
				IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(3,  9,  0,  1),   // j >= -9
                        										IV_(3, -9,  0, -1)), 1,
						new IHierarLeafPart(s)
				)
			)
		));

		// Assert we have the correct full iteration spaces
		/*
		ASSERT(s.getFullIterationSpace() == IPD(1, IP_(4, IV_(3,  9,  1,  0),
												           IV_(3, -9, -1,  0),
												           IV_(3,  9,  0,  1),
												           IV_(3, -9,  0, -1))));

		ASSERT(s.getFullIterationSpace() != IPD(1, IP_(4, IV_(3, 10,  1,  0),
												           IV_(3, -9, -1,  0),
												           IV_(3,  9,  0,  1),
												           IV_(3, -9,  0, -1))));

		ASSERT(s.getFullIterationSpace() != IPD(1, IP_(4, IV_(3,  9,  1,  0),
												           IV_(3, -9, -1,  0),
												           IV_(3,  9,  0,  1),
												           IV_(3, -8,  0, -1))));
		*/

		// Perform partitioning
		DebugOutStream ofs =  DebugOutStream("Yu_1");
		comp.osHtml = &ofs;

		comp.computeImplementation();
		//comp.performAffineSetSpacePartitioning();
		//comp.performAffineLatticeSpacePartitioning();
	}

	/* LU - Bondhugula's benchmark
	 * for (k = 0; k < N; k++) {
     *   for (j = k + 1; j < N; j++) {
     *     a[k][j] = a[k][j]/a[k][k];
     *   }
     *   for (i = k + 1; i < N; i++) {
     *     for (j = k + 1; j < N; j++)   {
     *       a[i][j] = a[i][j] - a[i][k]*a[k][j];
     *     }
     *   }
     * }
	 */
	void LU_Bondhugula() {
		IComputation comp = setupComputation(1);

		IStatement& sA = comp.addStatement();
		IStatement& sB = comp.addStatement();

		comp.setOriginalExecutionOrder(IHierarPart::fromLexicoLinearOrderedSubNodes(comp, IPDC(false, 1, IV_(2,  0,  1)), 1,
			IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(3,  0,  0,  1),
                    										IV_(3, -1,  1, -1)), 2,
				IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(4, -1,  0, -1,  1),
                         										IV_(4, -1,  1,  0, -1)), 1,
					new IHierarLeafPart(sA)
				),
				IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(4, -1,  0, -1,  1),
																IV_(4, -1,  1,  0, -1)), 1,
					IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(5, -1,  0, -1,  0,  1),
                             										IV_(5, -1,  1,  0,  0, -1)), 1,
						new IHierarLeafPart(sB)
					)
				)
			)
		));

		comp.addReference(sA, *v[0], WRITE, IAT(2, IV_(4,  0,  0,  1,  0),   // k
			       								   IV_(4,  0,  0,  0,  1))); // j

		comp.addReference(sA, *v[0], READ,  IAT(2, IV_(4,  0,  0,  1,  0),   // k
                							       IV_(4,  0,  0,  0,  1))); // j

		comp.addReference(sA, *v[0], READ,  IAT(2, IV_(4,  0,  0,  1,  0),   // k
			       								   IV_(4,  0,  0,  1,  0))); // k

		comp.addReference(sB, *v[0], WRITE, IAT(2, IV_(5,  0,  0,  0,  1,  0),   // i
				                                   IV_(5,  0,  0,  0,  0,  1))); // j

		comp.addReference(sB, *v[0], READ,  IAT(2, IV_(5,  0,  0,  0,  1,  0),   // i
				                                   IV_(5,  0,  0,  0,  0,  1))); // j

		comp.addReference(sB, *v[0], READ,  IAT(2, IV_(5,  0,  0,  0,  1,  0),   // i
				                                   IV_(5,  0,  0,  1,  0,  0))); // k

		comp.addReference(sB, *v[0], READ,  IAT(2, IV_(5,  0,  0,  1,  0,  0),   // k
				                                   IV_(5,  0,  0,  0,  0,  1))); // j

		// Perform partitioning
		DebugOutStream ofs = DebugOutStream("LU_Bondhugula");
		comp.osHtml = &ofs;

		comp.performAffineSetSpacePartitioning();
		//comp.performAffineLatticeSpacePartitioning();


	}

	/* Cholesky - Ahmed
	 * for (k = 0; k < N; k++) {
	 *   a[k][k] = sqrt(a[k][k])
     *   for (i = k + 1; i < N; i++) {
     *     a[i][k] = a[i][k]/a[k][k];
     *     for (j = k + 1; j <= i; j++)   {
     *       a[i][j] = a[i][j] - a[i][k]*a[j][k];
     *     }
     *   }
     * }
	 */
	void Cholesky_Ahmed() {
		IComputation comp = setupComputation(1);

		IStatement& sqrt = comp.addStatement();
		IStatement& div = comp.addStatement();
		IStatement& macc = comp.addStatement();

		comp.setOriginalExecutionOrder(IHierarPart::fromLexicoLinearOrderedSubNodes(comp, IPDC(false, 1, IV_(2,  0,  1)), 1,
			IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(3,  0,  0,  1),
                    										IV_(3, -1,  1, -1)), 2,
				new IHierarLeafPart(sqrt),
				IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(4, -1,  0, -1,  1),
                         										IV_(4, -1,  1,  0, -1)), 2,
					new IHierarLeafPart(div),
					IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(5, -1,  0, -1,  0,  1),
                             										IV_(5,  0,  0,  0,  1, -1)), 1,
                        new IHierarLeafPart(macc)
					)
				)
			)
		));

		comp.addReference(sqrt, *v[0], WRITE, IAT(2, IV_(3,  0,  0,  1),   // k
			       								     IV_(3,  0,  0,  1))); // k

		comp.addReference(sqrt, *v[0], READ, IAT(2, IV_(3,  0,  0,  1),   // k
			       								    IV_(3,  0,  0,  1))); // k


		comp.addReference(div, *v[0], WRITE, IAT(2, IV_(4,  0,  0,  0,  1),   // i
                							        IV_(4,  0,  0,  1,  0))); // k

		comp.addReference(div, *v[0], READ,  IAT(2, IV_(4,  0,  0,  0,  1),   // i
                							        IV_(4,  0,  0,  1,  0))); // k

		comp.addReference(div, *v[0], READ,  IAT(2, IV_(4,  0,  0,  1,  0),   // k
                							        IV_(4,  0,  0,  1,  0))); // k


		comp.addReference(macc, *v[0], WRITE, IAT(2, IV_(5,  0,  0,  0,  1,  0),   // i
			       							  	     IV_(5,  0,  0,  0,  0,  1))); // j

		comp.addReference(macc, *v[0], READ,  IAT(2, IV_(5,  0,  0,  0,  1,  0),   // i
			       							  	     IV_(5,  0,  0,  0,  0,  1))); // j

		comp.addReference(macc, *v[0], READ,  IAT(2, IV_(5,  0,  0,  0,  1,  0),   // i
			       							  	     IV_(5,  0,  0,  1,  0,  0))); // k

		comp.addReference(macc, *v[0], READ,  IAT(2, IV_(5,  0,  0,  0,  1,  0),   // j
			       							  	     IV_(5,  0,  0,  1,  0,  0))); // k



		// Perform partitioning
		DebugOutStream ofs = DebugOutStream("Cholesky_Ahmed");
		comp.osHtml = &ofs;

		comp.performAffineSetSpacePartitioning();
		//comp.performAffineLatticeSpacePartitioning();


	}

	/* Bifurcation
	 * for (i = 1; i < N; i++) {
     *   for (j = 0; j < N - i; j++) {
     *     for (k = 0; k < i; k++)   {
     * 		 v[j][i] = v[j][i] + v[j][k] * v[j + k + 1][i - k - 1]
     * } } }
	 */
	void bifurcation_Infernal() {
		IComputation comp = setupComputation(1);

		IStatement& s = comp.addStatement();

		comp.setOriginalExecutionOrder(IHierarPart::fromLexicoLinearOrderedSubNodes(comp, IPDC(false, 1, IV_(2,  0,  1)), 1,
			IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(3, -1,  0,  1),   // i >= 1
                    										            IV_(3, -1,  1, -1)), 1,
				IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(4,  0,  0,  0,  1),   // j >= 0
                        										            IV_(4, -1,  1, -1, -1)), 1,
					IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(5,  0,  0,  0,  0,  1), // k >= 0
                            										            IV_(5, -1,  0,  1,  0, -1)), 1,
                        new IHierarLeafPart(s)
					)
				)
			)
		));

		comp.addReference(s, *v[0], WRITE, IAT(2, IV_(5,  0,  0,  0,  1,  0),   // j
			       							  	  IV_(5,  0,  0,  1,  0,  0))); // i

		comp.addReference(s, *v[0], READ,  IAT(2, IV_(5,  0,  0,  0,  1,  0),   // j
			       							  	  IV_(5,  0,  0,  1,  0,  0))); // i

		comp.addReference(s, *v[0], READ,  IAT(2, IV_(5,  0,  0,  0,  1,  0),   // j
			       							  	  IV_(5,  0,  0,  0,  0,  1))); // k

		comp.addReference(s, *v[0], READ,  IAT(2, IV_(5,  1,  0,  0,  1,  1),   // j + k + 1
			       							  	  IV_(5, -1,  0,  1,  0, -1))); // i - k - 1

		// Assert we have the correct full iteration spaces
		/* Skipped: Equality testing not yet functional.
		*/

		// Perform partitioning
		DebugOutStream ofs = DebugOutStream("Bifurcation_Infernal");
		comp.osHtml = &ofs;

		comp.computeImplementation();
		ASSERT_EQ(*comp.spaceSolutions, IH(G, true, 2, IV_(5, 1, 0, 0, 0, 0),
				                                       IV_(5, 0, 1, 0, 0, 0)));

		ASSERT_EQ(*comp.timeSolutions, IH(G, true, 2, IV_(5, 0, 1, 0, 0, 0),
                                                      IV_(5, 1, 0, 0, 0, 0))
                                       + IH(G, false, 2, IV_(5, 0, 0, 1, 1, 0),
                                    		             IV_(5, 0, 0, 0, -1, 0)));

		//comp.performAffineSetSpacePartitioning();
		//comp.performAffineLatticeSpacePartitioning();

		//comp.parallellize<ICone>();
		//comp.performAffineSetSpacePartitioning();
		////comp.performAffineLatticeSpacePartitioning();


	}

	/* Modified bifurcation
	 * for (i = 1; i < N; i++) {
     *   for (k = 0; k < i; k++)   {
     * 	   v[i] = v[i] ++ v[k]
     * } }
	 */
	void A() {
		IComputation comp = setupComputation(1);

		IStatement& s = comp.addStatement();

		comp.setOriginalExecutionOrder(IHierarPart::fromLexicoLinearOrderedSubNodes(comp, IPDC(false, 1, IV_(2,  0,  1)), 1,
			IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(3, -1,  0,  1),
                    										IV_(3, -1,  1, -1)), 1,
				IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(4,  0,  0,  0,  1),
                        										IV_(4, -1,  0,  1, -1)), 1,
                    new IHierarLeafPart(s)
				)
			)
		));

		comp.addReference(s, *v[0], WRITE, IAT(1, IV_(4,  0,  0,  1,  0))); // i
		comp.addReference(s, *v[0], READ,  IAT(1, IV_(4,  0,  0,  1,  0))); // i
		comp.addReference(s, *v[0], READ,  IAT(1, IV_(4,  0,  0,  0,  1))); // k

		// Assert we have the correct full iteration spaces
		/* Skipped: Equality testing not yet functional.
		*/

		// Perform partitioning
		DebugOutStream ofs = DebugOutStream("A");
		comp.osHtml = &ofs;

		comp.computeImplementation();

		ASSERT_EQ(*comp.spaceSolutions, IH(G, true, 2, IV_(4, 1, 0, 0, 0),
				                                       IV_(4, 0, 1, 0, 0)));

		ASSERT_EQ(*comp.timeSolutions, IH(G, true, 2, IV_(4, 0, 1, 0, 0),
                                                      IV_(4, 1, 0, 0, 0))
                                       + IH(G, false, 1, IV_(4, 0, 0, 1, 0)));
		//comp.performAffineSetSpacePartitioning();
		//comp.performAffineLatticeSpacePartitioning();

		//comp.parallellize<ICone>();
		//comp.performAffineSetSpacePartitioning();
		////comp.performAffineLatticeSpacePartitioning();


	}

	/* Modified bifurcation
	 * for (i = 1; i < N; i++) {
     *   for (j = 0; j < N - i; j++) {
     *     for (k = 0; k < i; k++)   {
     * 		 v[j][i] = v[j][i] + v[j][k] * v[j + i - k][k]
     * } } }
	 */
	void bifurcation_Modified() {
		IComputation comp = setupComputation(1);

		IStatement& s = comp.addStatement();

		comp.setOriginalExecutionOrder(IHierarPart::fromLexicoLinearOrderedSubNodes(comp, IPDC(false, 1, IV_(2,  0,  1)), 1,
			IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(3, -1,  0,  1),   // i >= 1
                    										            IV_(3, -1,  1, -1)), 1,
				IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(4,  0,  0,  0,  1),   // j >= 0
                        										            IV_(4, -1,  1, -1, -1)), 1,
					IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(5,  0,  0,  0,  0,  1), // k >= 0
                            										            IV_(5, -1,  0,  1,  0, -1)), 1,
                        new IHierarLeafPart(s)
					)
				)
			)
		));

		comp.addReference(s, *v[0], WRITE, IAT(2, IV_(5,  0,  0,  0,  1,  0),   // j
			       							  	     IV_(5,  0,  0,  1,  0,  0))); // i

		comp.addReference(s, *v[0], READ,  IAT(2, IV_(5,  0,  0,  0,  1,  0),   // j
			       							  	     IV_(5,  0,  0,  1,  0,  0))); // i

		comp.addReference(s, *v[0], READ,  IAT(2, IV_(5,  0,  0,  0,  1,  0),   // j
			       							  	     IV_(5,  0,  0,  0,  0,  1))); // k

		comp.addReference(s, *v[0], READ,  IAT(2, IV_(5,  0,  0,  1,  1, -1),   // j + i - k
			       							  	     IV_(5,  0,  0,  0,  0,  1))); // k

		// Assert we have the correct full iteration spaces
		/* Skipped: Equality testing not yet functional.
		*/

		// Perform partitioning
		DebugOutStream ofs = DebugOutStream("Bifurcation_Modified");
		comp.osHtml = &ofs;

		comp.computeImplementation();

		ASSERT_EQ(*comp.spaceSolutions, IH(G, true, 2, IV_(5, 1, 0, 0, 0, 0),
				                                       IV_(5, 0, 1, 0, 0, 0)));

		ASSERT_EQ(*comp.timeSolutions, IH(G, true, 2, IV_(5, 1, 0, 0, 0, 0),
                                                      IV_(5, 0, 1, 0, 0, 0))
                                       + IH(G, false, 2, IV_(5, 0, 0, 1, 1, 0),
                                                         IV_(5, 0, 0, 0, -1, 0)));


	}

	/* Pair bifurcation (consan)
	 * for (dx = 0; dx <= X; dx++) {
	 *	   for (j = dx - 1; j <= X - 1; j++) {
	 *		   i = j - dx + 1;
	 *         for (dy = 0; dy <= Y; dy++) {
	 *      	   for (l = dy - 1; l <= Y - 1; l++) {
     * 			       k = l - dy + 1;
	 *      	       for (a = j - dx + 2; a < j - 1; a++) {
	 *      			   for (b = l - dy + 2; b < l - 1; b++) {
	 *   	                   v[j][dx][l][dy] = v[j][dx][l][dy] +
	 *                             v[j - 1][j - a - 1][l - 1][l - b - 1]
	 * 		                       * v[a - 1][a - j + dx - 1][b - 1][b - l + dy - 1]
	 * }   }   }   }   }   }
	 */
	void pairBifurcation_Consan() {
		IComputation comp = setupComputation(1);

		IStatement& s = comp.addStatement();

		comp.setOriginalExecutionOrder(IHierarPart::fromLexicoLinearOrderedSubNodes(comp,
			IPDC(false, 2, IV_(3,  0,  1,  0),   // X >= 0
                     IV_(3,  0,  0,  1)), 1,
			IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(4,  0,  0,  0,  1),
                    										IV_(4,  0,  1,  0, -1)), 1,
				IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(5,  1,  0,  0, -1,  1),
                        										IV_(5, -1,  1,  0,  0, -1)), 1,
					IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(6,  0,  0,  0,  0,  0,  1),
                             										IV_(6,  0,  0,  1,  0,  0, -1)), 1,
						IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(7,  1,  0,  0,  0,  0, -1,  1),
                                										IV_(7, -1,  0,  1,  0,  0,  0, -1)), 1,
							IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(8, -2,  0,  0,  1, -1,  0,  0,  1),
                                    										IV_(8, -2,  0,  0,  0,  1,  0,  0, -1)), 1,
								IHierarPart::fromLexicoLinearOrderedSubNodes(IPDC(false, 2, IV_(9, -2,  0,  0,  0,  0,  1, -1,  0,  1),
                                        										IV_(9, -2,  0,  0,  0,  0,  0,  1,  0, -1)), 1,
			                        new IHierarLeafPart(s)
								)
							)
						)
					)
				)
			)
		));

		comp.addReference(s, *v[0], WRITE, IAT(4, IV_(9,  0,  0,  0,  0,  1,  0,  0,  0,  0),   // j
			       								     IV_(9,  0,  0,  0,  1,  0,  0,  0,  0,  0),   // dx
			       								     IV_(9,  0,  0,  0,  0,  0,  0,  1,  0,  0),   // l
			       								     IV_(9,  0,  0,  0,  0,  0,  1,  0,  0,  0))); // dy

		comp.addReference(s, *v[0], READ,  IAT(4, IV_(9,  0,  0,  0,  0,  1,  0,  0,  0,  0),   // j
			       								     IV_(9,  0,  0,  0,  1,  0,  0,  0,  0,  0),   // dx
			       								     IV_(9,  0,  0,  0,  0,  0,  0,  1,  0,  0),   // l
			       								     IV_(9,  0,  0,  0,  0,  0,  1,  0,  0,  0))); // dy

		comp.addReference(s, *v[0], READ,  IAT(4, IV_(9, -1,  0,  0,  0,  1,  0,  0,  0,  0),   // j - 1
			       								     IV_(9, -1,  0,  0,  0,  1,  0,  0, -1,  0),   // j - a - 1
			       								     IV_(9, -1,  0,  0,  0,  0,  0,  1,  0,  0),   // l - 1
			       								     IV_(9, -1,  0,  0,  0,  0,  0,  1,  0, -1))); // l - b - 1

		comp.addReference(s, *v[0], READ,  IAT(4, IV_(9, -1,  0,  0,  0,  0,  0,  0,  1,  0),   // a - 1
			       								     IV_(9, -1,  0,  0,  1, -1,  0,  0,  1,  0),   // a - j + dx - 1
			       								     IV_(9, -1,  0,  0,  0,  0,  0,  0,  0,  1),   // b - 1
			       								     IV_(9, -1,  0,  0,  0,  0,  1, -1,  0,  1))); // b - l + dy - 1

		// Assert we have the correct full iteration spaces
		/* Skipped: Equality testing not yet functional.
		*/

		// Perform partitioning
		DebugOutStream ofs = DebugOutStream("Pairwise bifurcation (Consan)");
		comp.osHtml = &ofs;

		comp.computeImplementation();
		//comp.performAffineSetSpacePartitioning();
		//comp.performAffineLatticeSpacePartitioning();
		//comp.parallellize<ICone>();
		//comp.performAffineSetSpacePartitioning();
		////comp.performAffineLatticeSpacePartitioning();


	}



	/*void CComputation_SpacePartition_2DUniformLattice_YuRot() {
		IComputation comp = setupComputation(2, 2, 0);

		s[0]->setIterationSpace(IP(4, IV_(3,  1,  0,  9),
							 	         IV_(3, -1,  0, -9),
							 	         IV_(3,  0,  1,  8),
							 	         IV_(3,  0, -1, -8)));
		s[1]->setIterationSpace(IP(4, IV_(3,  1,  0,  9),
							 	   	IV_(3, -1,  0, -),
							 	   	IV_(3,  0,  1,  8),
							 	   IV_(3,  0, -1, -8)));

		comp.addReference(*s[0], *v[0], WRITE, IAT(2, IV_(3,  4, -1,  3),
				                                      IV_(3,  2,  1, -2)));
		comp.addReference(*s[0], *v[0], READ,  IAT(2, IV_(3,  1,  1, -1),
										              IV_(3,  1, -1,  2)));
		comp.addReference(*s[1], *v[0], READ,  IAT(2, IV_(3,  1,  0, 2),
				                                      IV_(3,  1,  0, -2)));
		comp.addReference(*s[1], *v[1], WRITE, IAT(2, IV_(3,  1,  0,  0),
				                                      IV_(3,  0,  1,  0)));
		comp.addReference(*s[1], *v[1], READ,  IAT(2, IV_(3,  1,  0,  0),
                				                      IV_(3,  0,  1, -2)));

		DebugOutStream ofs = DebugOutStream("CComputation_SpacePartition_2DUniformLattice_YuRot");
		comp.osHtml = &ofs;

		comp.performAffineSetSpacePartitioning();
		comp.performAffineLatticeSpacePartitioning();


	}*/
}

cute::suite* Test_Computation_NewHadron_runSuite(){
	cute::suite &s = *(new cute::suite("Computation_NewHadron"));

	s.push_back(CUTE(Test_Computation_NewHadron::LexiCone));

	s.push_back(CUTE(Test_Computation_NewHadron::Trivial1DParallelism_2S_1P_B));
	s.push_back(CUTE(Test_Computation_NewHadron::Trivial1DParallelism_2S_1P_C));

	s.push_back(CUTE(Test_Computation_NewHadron::Trivial1DSequential_While));
	s.push_back(CUTE(Test_Computation_NewHadron::Trivial1DSequential_While_L));
	s.push_back(CUTE(Test_Computation_NewHadron::Trivial1DSequential));

	s.push_back(CUTE(Test_Computation_NewHadron::simple2DUniform_Unbounded));
	s.push_back(CUTE(Test_Computation_NewHadron::simple2DUniform));

	s.push_back(CUTE(Test_Computation_NewHadron::Trivial1DParallelism));
	s.push_back(CUTE(Test_Computation_NewHadron::Trivial1DParallelism_2S));
	s.push_back(CUTE(Test_Computation_NewHadron::Trivial1DParallelism_2S_1P));
	s.push_back(CUTE(Test_Computation_NewHadron::Trivial1DParallelism_2S_1P_B));

	s.push_back(CUTE(Test_Computation_NewHadron::Trivial2DParallelism_SingleStatement));
	s.push_back(CUTE(Test_Computation_NewHadron::Trivial2DParallelism));
	s.push_back(CUTE(Test_Computation_NewHadron::Trivial2DParallelism_B));
	//s.push_back(CUTE(Test_Computation_NewHadron::Trivial2DParallelism_C));

	s.push_back(CUTE(Test_Computation_NewHadron::A));

	s.push_back(CUTE(Test_Computation_NewHadron::bifurcation_Modified));
	s.push_back(CUTE(Test_Computation_NewHadron::bifurcation_Infernal));

	s.push_back(CUTE(Test_Computation_NewHadron::pairBifurcation_Consan));

	s.push_back(CUTE(Test_Computation_NewHadron::LU_Bondhugula));
	s.push_back(CUTE(Test_Computation_NewHadron::Cholesky_Ahmed));

	s.push_back(CUTE(Test_Computation_NewHadron::Yu_1));


	//s.push_back(CUTE(Test_Computation::CComputation_SpacePartition_2DUniformLattice_YuRot));

	return &s;
}
