#ifdef OBSOLETE

#include "../cute/cute.h"
#include "../cute/cute_suite.h"

#include "partyql/graph/Graph.h"
#include "partyql/basicmath/scalar/Integer.h"
#include "partyql/basicmath/Vektor.h"
#include "partyql/basicmath/Module.h"
#include "partyql/basicmath/SetRelation.h"
#include "partyql/basicmath/AffineTransformation.h"
#include "partyql/basicmath/PolyhedralDomain.h"
#include "partyql/RDMGraph.h"
#include "partyql/HierarGraph.h"
#include "partyql/Computation.h"
#include "partyql/Statement.h"

#include "utils/SVG.h"
#include "utils/SVGELattice.h"
#include "utils/SVGViewbox.h"
#include "utils/ToString.h"
#include "utils/Html.h"

#include "partyql/SpacePartitioner.h"

#include "Test_Funcs.h"

#include <fstream>
using std::ofstream;

// #define ENABLE_PDF

using namespace AlgoTrans;

namespace Test_Computation {

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
	typedef CRDMGraph<I> RDMGraph;

	typedef CVertex<CPolyheder<I>, ILatticeRelation> RADLatticeVertex;
	typedef CEdge<CPolyheder<I>, ILatticeRelation> RADLatticeEdge;
	typedef CRADGraph<I, ILatticeRelation> RADLatticeGraph;

	typedef CVertex<CPolyheder<I>, IModuleRelation> RADModuleVertex;
	typedef CEdge<CPolyheder<I>, IModuleRelation> RADModuleEdge;
	typedef CRADGraph<I, IModuleRelation> RADModuleGraph;

	typedef CHierarNode<CComputation<I>, IPolyhedralDomain, CStatement<I> > IHierarNode;
	//typedef CHierarPart<CComputation<I>, IPolyhedralDomain, CStatement<I> > IHierarNode;
	//typedef CHierarLeafPart<CComputation<I>, IPolyhedralDomain, CStatement<I> > IHierarNode;

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

	ofstream* ofsAll;

	ofstream* initXHTMLStream(string fileName) {
		ofstream* ofs = new ofstream(("examples/" + fileName + ".xhtml").c_str());

		*ofs << "<html xmlns=\"http://www.w3.org/1999/xhtml\"><body>";
		//*ofsAll << "<a href=\"" + fileName + ".xhtml\">" + fileName + "</a><br>";

		return ofs;
	}

	void closeXHTMLStream(ofstream* ofs) {
		*ofs << "</body></html>";
		ofs->close();
	}

	/* 1D Sequential loop
	 * for ( i = 0 ;   ; i++) {
	 *     A[i] := A[i] + A[i - 1];
	 * }	 */
	void Trivial1DSequential_While() {
		IComputation comp = setupComputation(1);

		IStatement& s = comp.addStatement();

		comp.setOriginalExecutionOrder(IHierarNode::fromLexicoLinearOrderedSubNodes(comp, pIP_D(0, 0), 1,
				IHierarNode::fromLexicoLinearOrderedSubNodes(pIP_D(1, IV_(2,  0,  1)), 1,
						new IHierarNode(&s)
				)
		));

		IAffineTransformation iat_i = IAT(1, IV_(2,  0,  1));
		IAffineTransformation iat_i_1 = IAT(1, IV_(2, -1,  1));
		comp.addReference(s, *v[0], WRITE, iat_i);
		comp.addReference(s, *v[0], READ,  iat_i);
		comp.addReference(s, *v[0], READ,  iat_i_1);

		// Assert we have the correct full iteration spaces
		ASSERT(s.getFullIterationSpace() == IPD(1, IP_(1, IV_(2,  0,  1))));

		// Perform partitioning
		ofstream* ofs = initXHTMLStream("Trivial1DSequential_While");
		comp.osHtml = ofs;

		comp.computeImplementation<CFlat<IModule> >();

		//comp.performAffineSetSpacePartitioning();
		//comp.performAffineLatticeSpacePartitioning();

		closeXHTMLStream(ofs);

		ASSERT_EQ(comp.dependenceFlatClosure->getUniqueEdgeBetween(0, 0).getData(),
			      IMR_Con(1, 1, 1, 0));

		comp.dependenceFlatClosure->getUniqueEdgeBetween(0, 0).getData().print(true);
		//comp.dependenceLatticeClosure->getUniqueEdgeBetween(0, 0).getData().print(true);
	}

	/* 1D Sequential loop
	 * for ( i = 0 ;   ; i++) {
	 *     A[i] := A[i] + A[i - 2];
	 * }	 */
	void Trivial1DSequential_While_L() {
		IComputation comp = setupComputation(1);

		IStatement& s = comp.addStatement();

		comp.setOriginalExecutionOrder(IHierarNode::fromLexicoLinearOrderedSubNodes(comp, pIP_D(0, 0), 1,
				IHierarNode::fromLexicoLinearOrderedSubNodes(pIP_D(1, IV_(2,  0,  1)), 1,
						new IHierarNode(&s)
				)
		));

		IAffineTransformation iat_i = IAT(1, IV_(2,  0,  1));
		IAffineTransformation iat_i_2 = IAT(1, IV_(2, -2,  1));
		comp.addReference(s, *v[0], WRITE, iat_i);
		comp.addReference(s, *v[0], READ,  iat_i);
		comp.addReference(s, *v[0], READ,  iat_i_2);

		// Assert we have the correct full iteration spaces
		ASSERT(s.getFullIterationSpace() == IPD(1, IP_(1, IV_(2,  0,  1))));

		// Perform partitioning
		ofstream* ofs = initXHTMLStream("Trivial1DSequential_While_L");
		comp.osHtml = ofs;

		comp.performAffineSetSpacePartitioning();
		//comp.performAffineLatticeSpacePartitioning();

		closeXHTMLStream(ofs);

		ASSERT_EQ(comp.dependenceFlatClosure->getUniqueEdgeBetween(0, 0).getData(),
			      IMR_Con(1, 1, 1, 0));

		comp.dependenceFlatClosure->getUniqueEdgeBetween(0, 0).getData().print(true);
		//comp.dependenceLatticeClosure->getUniqueEdgeBetween(0, 0).getData().print(true);
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

		comp.setOriginalExecutionOrder(IHierarNode::fromLexicoLinearOrderedSubNodes(comp, pIP_D(1, IV_(2, -1,  1)), 1,
				IHierarNode::fromLexicoLinearOrderedSubNodes(pIP_D(2, IV_(3,  0,  0,  1),
                        										IV_(3,  0,  1, -1)), 1,
						new IHierarNode(&s)
				)
		));

		IAffineTransformation iat_i = IAT(1, IV_(3,  0,  0,  1));
		IAffineTransformation iat_i_1 = IAT(1, IV_(3, -1,  0,  1));
		comp.addReference(s, *v[0], WRITE, iat_i);
		comp.addReference(s, *v[0], READ,  iat_i);
		comp.addReference(s, *v[0], READ,  iat_i_1);

		// Assert we have the correct full iteration spaces
		ASSERT_EQ(s.getFullIterationSpace(), IPD(1, IP_(3, IV_(3, -1,  1,  0),
											 	            IV_(3,  0,  0,  1),
												            IV_(3,  0,  1, -1))));

		// Perform partitioning
		ofstream* ofs = initXHTMLStream("Trivial1DSequential");
		comp.osHtml = ofs;

		comp.performAffineSetSpacePartitioning();
		//comp.performAffineLatticeSpacePartitioning();

		closeXHTMLStream(ofs);

		ASSERT_EQ(comp.dependenceFlatClosure->getUniqueEdgeBetween(0, 0).getData(),
			      IMR_Con(1, 1, 2, 0));
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

		comp.setOriginalExecutionOrder(IHierarNode::fromLexicoLinearOrderedSubNodes(comp, pIP_D(0, 0), 1,
			IHierarNode::fromLexicoLinearOrderedSubNodes(pIP_D(2, IV_(2,  0,  1),
                      										IV_(2, 77, -1)), 1,
				IHierarNode::fromLexicoLinearOrderedSubNodes(pIP_D(2, IV_(3,  0,  0,  1),
                        										IV_(3, 77,  0, -1)), 1,
						new IHierarNode(&s)
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
		ofstream* ofs = initXHTMLStream("simple2DUniform");
		comp.osHtml = ofs;

		comp.performAffineSetSpacePartitioning();
		//comp.performAffineLatticeSpacePartitioning();

		closeXHTMLStream(ofs);

		ASSERT_EQ(comp.dependenceFlatClosure->getUniqueEdgeBetween(0, 0).getData(),
			      IMR_Con(2, 2, 1, 0));
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

		comp.setOriginalExecutionOrder(IHierarNode::fromLexicoLinearOrderedSubNodes(comp, pIP_D(0, 0), 1,
			IHierarNode::fromLexicoLinearOrderedSubNodes(pIP_D(1, IV_(2,  0,  1)), 1,
				IHierarNode::fromLexicoLinearOrderedSubNodes(pIP_D(1, IV_(3,  0,  0,  1)), 1,
						new IHierarNode(&s)
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
		ofstream* ofs = initXHTMLStream("simple2DUniform_Unbounded");
		comp.osHtml = ofs;

		comp.performAffineSetSpacePartitioning();
		//comp.performAffineLatticeSpacePartitioning();

		closeXHTMLStream(ofs);

		ASSERT_EQ(comp.dependenceFlatClosure->getUniqueEdgeBetween(0, 0).getData(),
			      IMR_Con(2, 2, 1, 0));
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
	void Trivial2DParallelism() {
		IComputation comp = setupComputation(2);

		IStatement& sC = comp.addStatement();
		IStatement& sD = comp.addStatement();

		comp.setOriginalExecutionOrder(IHierarNode::fromLexicoLinearOrderedSubNodes(comp,
			pIP_D(2, IV_(3, -1,  1,  0),   // n >= 1
                     IV_(3, -2,  0,  1)), 1,
			IHierarNode::fromLexicoLinearOrderedSubNodes(pIP_D(2, IV_(4,  0,  0,  0,  1),   // i >= 0
                    										IV_(4,  0,  1,  0, -1)), 1,
				IHierarNode::fromLexicoLinearOrderedSubNodes(pIP_D(2, IV_(5,  0,  0,  0,  0,  1),   // j >= 0
                        										IV_(5,  0,  0,  1,  0, -1)), 2,
						new IHierarNode(&sC),
						new IHierarNode(&sD)
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
		ASSERT_EQ(sC.getFullIterationSpace(), IPD(1, IP_(6, IV_(5, -1,  1,  0,  0,  0),
												           IV_(5, -2,  0,  1,  0,  0),
												           IV_(5,  0,  0,  0,  1,  0),
												           IV_(5,  0,  1,  0, -1,  0),
												           IV_(5,  0,  0,  0,  0,  1),
												           IV_(5,  0,  0,  1,  0, -1))));

		ASSERT_EQ(sD.getFullIterationSpace(), IPD(1, IP_(6, IV_(5, -1,  1,  0,  0,  0),
														   IV_(5, -2,  0,  1,  0,  0),
														   IV_(5,  0,  0,  0,  1,  0),
														   IV_(5,  0,  1,  0, -1,  0),
														   IV_(5,  0,  0,  0,  0,  1),
														   IV_(5,  0,  0,  1,  0, -1))));

		// Perform partitioning
		ofstream* ofs = initXHTMLStream("Trivial2DParallelism");
		comp.osHtml = ofs;

		comp.computeImplementation<CFlat<IModule > >();

		comp.performAffineSetSpacePartitioning();
		//comp.performAffineLatticeSpacePartitioning();

		closeXHTMLStream(ofs);

		/*
		ASSERT_EQ(comp.dependenceFlatClosure->getUniqueEdgeBetween(0, 0).getData(),
			      IMR_Con(2, 2, 3, 2, IV_(7,  0,  0,  0,  1,  0, -1,  0),
						              IV_(7,  0,  0,  0,  0,  1,  0, -1)));
   		ASSERT_EQ(comp.dependenceFlatClosure->getUniqueEdgeBetween(1, 1).getData(),
   			      IMR_Con(2, 2, 3, 2, IV_(7,  0,  0,  0,  1,  0, -1,  0),
   						              IV_(7,  0,  0,  0,  0,  1,  0, -1)));
   		ASSERT_EQ(comp.dependenceFlatClosure->getUniqueEdgeBetween(0, 1).getData(),
   			      IMR_Con(2, 2, 3, 2, IV_(7,  0,  0,  0,  1,  0, -1,  0),
   						              IV_(7,  0,  0,  0,  0,  1,  0, -1)));
   		ASSERT_EQ(comp.dependenceFlatClosure->getUniqueEdgeBetween(1, 0).getData(),
   			      IMR_Con(2, 2, 3, 2, IV_(7,  0,  0,  0,  1,  0, -1,  0),
  						              IV_(7,  0,  0,  0,  0,  1,  0, -1)));
  						              */
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

		comp.setOriginalExecutionOrder(IHierarNode::fromLexicoLinearOrderedSubNodes(comp, pIP_D(0, 0), 1,
			IHierarNode::fromLexicoLinearOrderedSubNodes(pIP_D(2, IV_(2,  9,  1),   // i >= -9
                    										IV_(2, -9, -1)), 1,
				IHierarNode::fromLexicoLinearOrderedSubNodes(pIP_D(2, IV_(3,  9,  0,  1),   // j >= -9
                        										IV_(3, -9,  0, -1)), 1,
						new IHierarNode(&s)
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
		ofstream* ofs = initXHTMLStream("Yu_1");
		comp.osHtml = ofs;

		comp.performAffineSetSpacePartitioning();
		//comp.performAffineLatticeSpacePartitioning();

		closeXHTMLStream(ofs);
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

		comp.setOriginalExecutionOrder(IHierarNode::fromLexicoLinearOrderedSubNodes(comp, pIP_D(1, IV_(2,  0,  1)), 1,
			IHierarNode::fromLexicoLinearOrderedSubNodes(pIP_D(2, IV_(3,  0,  0,  1),
                    										IV_(3, -1,  1, -1)), 2,
				IHierarNode::fromLexicoLinearOrderedSubNodes(pIP_D(2, IV_(4, -1,  0, -1,  1),
                         										IV_(4, -1,  1,  0, -1)), 1,
					new IHierarNode(&sA)
				),
				IHierarNode::fromLexicoLinearOrderedSubNodes(pIP_D(2, IV_(4, -1,  0, -1,  1),
																IV_(4, -1,  1,  0, -1)), 1,
					IHierarNode::fromLexicoLinearOrderedSubNodes(pIP_D(2, IV_(5, -1,  0, -1,  0,  1),
                             										IV_(5, -1,  1,  0,  0, -1)), 1,
						new IHierarNode(&sB)
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
		ofstream* ofs = initXHTMLStream("LU_Bondhugula");
		comp.osHtml = ofs;

		comp.performAffineSetSpacePartitioning();
		//comp.performAffineLatticeSpacePartitioning();

		closeXHTMLStream(ofs);
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

		comp.setOriginalExecutionOrder(IHierarNode::fromLexicoLinearOrderedSubNodes(comp, pIP_D(1, IV_(2,  0,  1)), 1,
			IHierarNode::fromLexicoLinearOrderedSubNodes(pIP_D(2, IV_(3,  0,  0,  1),
                    										IV_(3, -1,  1, -1)), 2,
				new IHierarNode(&sqrt),
				IHierarNode::fromLexicoLinearOrderedSubNodes(pIP_D(2, IV_(4, -1,  0, -1,  1),
                         										IV_(4, -1,  1,  0, -1)), 2,
					new IHierarNode(&div),
					IHierarNode::fromLexicoLinearOrderedSubNodes(pIP_D(2, IV_(5, -1,  0, -1,  0,  1),
                             										IV_(5,  0,  0,  0,  1, -1)), 1,
                        new IHierarNode(&macc)
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
		ofstream* ofs = initXHTMLStream("Cholesky_Ahmed");
		comp.osHtml = ofs;

		comp.performAffineSetSpacePartitioning();
		//comp.performAffineLatticeSpacePartitioning();

		closeXHTMLStream(ofs);
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

		comp.setOriginalExecutionOrder(IHierarNode::fromLexicoLinearOrderedSubNodes(comp, pIP_D(1, IV_(2,  0,  1)), 1,
			IHierarNode::fromLexicoLinearOrderedSubNodes(pIP_D(2, IV_(3, -1,  0,  1),   // i >= 1
                    										IV_(3, -1,  1, -1)), 1,
				IHierarNode::fromLexicoLinearOrderedSubNodes(pIP_D(2, IV_(4,  0,  0,  0,  1),   // j >= 0
                        										IV_(4, -1,  1, -1, -1)), 1,
					IHierarNode::fromLexicoLinearOrderedSubNodes(pIP_D(2, IV_(5,  0,  0,  0,  0,  1), // k >= 0
                            										IV_(5, -1,  0,  1,  0, -1)), 1,
                        new IHierarNode(&s)
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
		ofstream* ofs = initXHTMLStream("Bifurcation_Infernal");
		comp.osHtml = ofs;

		comp.performAffineSetSpacePartitioning();
		//comp.performAffineLatticeSpacePartitioning();

		//comp.parallellize<ICone>();
		//comp.performAffineSetSpacePartitioning();
		////comp.performAffineLatticeSpacePartitioning();

		closeXHTMLStream(ofs);
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

		comp.setOriginalExecutionOrder(IHierarNode::fromLexicoLinearOrderedSubNodes(comp, pIP_D(1, IV_(2,  0,  1)), 1,
			IHierarNode::fromLexicoLinearOrderedSubNodes(pIP_D(2, IV_(3, -1,  0,  1),
                    										IV_(3, -1,  1, -1)), 1,
				IHierarNode::fromLexicoLinearOrderedSubNodes(pIP_D(2, IV_(4,  0,  0,  0,  1),
                        										IV_(4, -1,  0,  1, -1)), 1,
                    new IHierarNode(&s)
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
		ofstream* ofs = initXHTMLStream("A");
		comp.osHtml = ofs;

		comp.performAffineSetSpacePartitioning();
		//comp.performAffineLatticeSpacePartitioning();

		//comp.parallellize<ICone>();
		//comp.performAffineSetSpacePartitioning();
		////comp.performAffineLatticeSpacePartitioning();

		closeXHTMLStream(ofs);
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

		comp.setOriginalExecutionOrder(IHierarNode::fromLexicoLinearOrderedSubNodes(comp, pIP_D(1, IV_(2,  0,  1)), 1,
			IHierarNode::fromLexicoLinearOrderedSubNodes(pIP_D(2, IV_(3, -1,  0,  1),   // i >= 1
                    										IV_(3, -1,  1, -1)), 1,
				IHierarNode::fromLexicoLinearOrderedSubNodes(pIP_D(2, IV_(4,  0,  0,  0,  1),   // j >= 0
                        										IV_(4, -1,  1, -1, -1)), 1,
					IHierarNode::fromLexicoLinearOrderedSubNodes(pIP_D(2, IV_(5,  0,  0,  0,  0,  1), // k >= 0
                            										IV_(5, -1,  0,  1,  0, -1)), 1,
                        new IHierarNode(&s)
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
		ofstream* ofs = initXHTMLStream("Bifurcation_Modified");
		comp.osHtml = ofs;

		comp.performAffineSetSpacePartitioning();
		//comp.performAffineLatticeSpacePartitioning();

		//comp.parallellize<ICone>();
		//comp.performAffineSetSpacePartitioning();
		////comp.performAffineLatticeSpacePartitioning();

		closeXHTMLStream(ofs);
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

		comp.setOriginalExecutionOrder(IHierarNode::fromLexicoLinearOrderedSubNodes(comp,
			pIP_D(2, IV_(3,  0,  1,  0),   // X >= 0
                     IV_(3,  0,  0,  1)), 1,
			IHierarNode::fromLexicoLinearOrderedSubNodes(pIP_D(2, IV_(4,  0,  0,  0,  1),
                    										IV_(4,  0,  1,  0, -1)), 1,
				IHierarNode::fromLexicoLinearOrderedSubNodes(pIP_D(2, IV_(5,  1,  0,  0, -1,  1),
                        										IV_(5, -1,  1,  0,  0, -1)), 1,
					IHierarNode::fromLexicoLinearOrderedSubNodes(pIP_D(2, IV_(6,  0,  0,  0,  0,  0,  1),
                             										IV_(6,  0,  0,  1,  0,  0, -1)), 1,
						IHierarNode::fromLexicoLinearOrderedSubNodes(pIP_D(2, IV_(7,  1,  0,  0,  0,  0, -1,  1),
                                										IV_(7, -1,  0,  1,  0,  0,  0, -1)), 1,
							IHierarNode::fromLexicoLinearOrderedSubNodes(pIP_D(2, IV_(8, -2,  0,  0,  1, -1,  0,  0,  1),
                                    										IV_(8, -2,  0,  0,  0,  1,  0,  0, -1)), 1,
								IHierarNode::fromLexicoLinearOrderedSubNodes(pIP_D(2, IV_(9, -2,  0,  0,  0,  0,  1, -1,  0,  1),
                                        										IV_(9, -2,  0,  0,  0,  0,  0,  1,  0, -1)), 1,
			                        new IHierarNode(&s)
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
		ofstream* ofs = initXHTMLStream("Pairwise bifurcation (Consan)");
		comp.osHtml = ofs;

		comp.performAffineSetSpacePartitioning();
		//comp.performAffineLatticeSpacePartitioning();
		//comp.parallellize<ICone>();
		//comp.performAffineSetSpacePartitioning();
		////comp.performAffineLatticeSpacePartitioning();

		closeXHTMLStream(ofs);
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

		ofstream* ofs = initXHTMLStream("CComputation_SpacePartition_2DUniformLattice_YuRot");
		comp.osHtml = ofs;

		comp.performAffineSetSpacePartitioning();
		comp.performAffineLatticeSpacePartitioning();

		closeXHTMLStream(ofs);
	}*/
}

cute::suite* Test_Computation_runSuite(){
	cute::suite &s = *(new cute::suite("Computation"));

	s.push_back(CUTE(Test_Computation::Trivial1DSequential_While));
	s.push_back(CUTE(Test_Computation::Trivial1DSequential_While_L));
	s.push_back(CUTE(Test_Computation::Trivial1DSequential));

	s.push_back(CUTE(Test_Computation::simple2DUniform));
	s.push_back(CUTE(Test_Computation::simple2DUniform_Unbounded));

	s.push_back(CUTE(Test_Computation::A));

	//s.push_back(CUTE(Test_Computation::bifurcation_Modified));
	//s.push_back(CUTE(Test_Computation::bifurcation_Infernal));

	//s.push_back(CUTE(Test_Computation::pairBifurcation_Consan));


	s.push_back(CUTE(Test_Computation::Trivial2DParallelism));
	s.push_back(CUTE(Test_Computation::Yu_1));
	s.push_back(CUTE(Test_Computation::LU_Bondhugula));
	s.push_back(CUTE(Test_Computation::Cholesky_Ahmed));


	//s.push_back(CUTE(Test_Computation::CComputation_SpacePartition_2DUniformLattice_YuRot));

	return &s;
}

#endif
