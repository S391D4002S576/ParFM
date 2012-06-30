#ifdef obsolete

/* OBSOLETE */

#include "cute.h"
#include "cute_suite.h"

#include "../partyql/graph/Graph.h"
#include "../partyql/basicmath/scalar/Integer.h"
#include "../partyql/basicmath/Vektor.h"
#include "../partyql/basicmath/Module.h"
#include "../partyql/basicmath/SetRelation.h"
#include "../partyql/basicmath/AffineTransformation.h"
#include "../partyql/RDMGraph.h"
#include "../partyql/Computation.h"

#include "../utils/SVG.h"
#include "../utils/SVGELattice.h"
#include "../utils/SVGViewbox.h"
#include "../utils/ToString.h"

#include "../partyql/SpacePartitioner.h"

#include "Test_Funcs.h"

#include <fstream>
using std::ofstream;

// #define ENABLE_PDF

using namespace AlgoTrans;

namespace Test_RDMGraph {
	typedef CInteger I;
	typedef CVector<I> IVector;
	typedef CAffineTransformation<I> IAffineTransformation;
	typedef CModuleRelation<I>::T IModuleRelation;
	typedef CLatticeRelation<I>::T ILatticeRelation;

	typedef CComputation<I> IComputation;
	typedef CVariable<I> IVariable;
	typedef CStatement<I> IStatement;

	int spaceFact = 25;

	I aiv1_0_0_0_0[]  = {I(1), I(0), I(0), I(0), I( 0)};
	I aiv0_1_0_0_0[]  = {I(0), I(1), I(0), I(0), I( 0)};
	I aiv0_0_0_0_1[]  = {I(0), I(0), I(0), I(0), I( 1)};
	I aiv1_0_0_0__1[] = {I(1), I(0), I(0), I(0), I(-1)};
	I aiv0_1_0_0__1[] = {I(0), I(1), I(0), I(0), I(-1)};
	IVector iv1_0_0_0_0, iv0_1_0_0_0, iv0_0_0_0_1, iv1_0_0_0__1, iv0_1_0_0__1;

	I aiv1__1_1__1_0_0__1[] = {I(1), I(-1), I(1), I(-1), I(0), I(0), I(-1)};
	I aiv1__1_1__1_0_0_1[]  = {I(1), I(-1), I(1), I(-1), I(0), I(0), I( 1)};
	I aiv1__1_1__1_0_0_0[]  = {I(1), I(-1), I(1), I(-1), I(0), I(0), I( 0)};
	I aiv1_0_1_0_0_0_0[]    = {I(1), I( 0), I(1), I( 0), I(0), I(0), I( 0)};
	I aiv0_1_0_1_0_0_0[]    = {I(0), I( 1), I(0), I( 1), I(0), I(0), I( 0)};
	IVector iv1__1_1__1_0_0__1, iv1__1_1__1_0_0_1, iv1__1_1__1_0_0_0, iv1_0_1_0_0_0_0, iv0_1_0_1_0_0_0;

	I aiv0_1_0_0_0_0__1[]    = {I( 0), I( 1), I( 0), I( 0), I( 0), I( 0), I(-1)};
	I aiv0_0_0_1_0_0_1[]     = {I( 0), I( 0), I( 0), I( 1), I( 0), I( 0), I( 1)};
	IVector iv0_1_0_0_0_0__1, iv0_0_0_1_0_0_1;

	I aiv0_1_1_0_0_0_0[]     = {I( 0), I( 1), I( 1), I( 0), I( 0), I( 0), I( 0)};
	I aiv1_0_0_0_0_0__1[]    = {I( 1), I( 0), I( 0), I( 0), I( 0), I( 0), I(-1)};
	IVector iv0_1_1_0_0_0_0, iv1_0_0_0_0_0__1;

	I aiv1_0_0_1_0_0_0[]     = {I( 1), I( 0), I( 0), I( 1), I( 0), I( 0), I( 0)};
	I aiv0_0_1_0_0_0_1[]     = {I( 0), I( 0), I( 1), I( 0), I( 0), I( 0), I( 1)};
	IVector iv1_0_0_1_0_0_0, iv0_0_1_0_0_0_1;

	I aiv1_0_0__6[]       = {I( 1), I( 0), I( 0), I(-6)};
	I aiv1_0_0__5[]       = {I( 1), I( 0), I( 0), I(-5)};
	I aiv1_0_0__1[]       = {I( 1), I( 0), I( 0), I(-1)};
	I aiv1_0_0_0[]        = {I( 1), I( 0), I( 0), I( 0)};
	IVector iv1_0_0__5, iv1_0_0__6, iv1_0_0__1, iv1_0_0_0;

	I aiv1_0_0_0__5[]        = {I( 1), I( 0), I( 0), I( 0), I(-5)};
	I aiv1_0_0_0__3[]        = {I( 1), I( 0), I( 0), I( 0), I(-3)};
	I aiv0_1_0_0__4[]        = {I( 0), I( 1), I( 0), I( 0), I(-4)};
	IVector iv1_0_0_0__5, iv1_0_0_0__3, iv0_1_0_0__4;

	void test_RDMGraph_initializeIVectors() {
		iv1_0_0_0_0 = IVector(5, aiv1_0_0_0_0);
		iv0_1_0_0_0 = IVector(5, aiv0_1_0_0_0);
		iv0_0_0_0_1 = IVector(5, aiv0_0_0_0_1);
		iv1_0_0_0__1 = IVector(5, aiv1_0_0_0__1);
		iv0_1_0_0__1 = IVector(5, aiv0_1_0_0__1);

		iv1_0_0_0 = IVector(4, aiv1_0_0_0);
		iv1_0_0__6 = IVector(4, aiv1_0_0__6);
		iv1_0_0__5 = IVector(4, aiv1_0_0__5);
		iv1_0_0__1 = IVector(4, aiv1_0_0__1);

		iv1_0_0_0__5 = IVector(5, aiv1_0_0_0__5);
		iv1_0_0_0__3 = IVector(5, aiv1_0_0_0__3);
		iv0_1_0_0__4 = IVector(5, aiv0_1_0_0__4);

		iv1__1_1__1_0_0__1 = IVector(7, aiv1__1_1__1_0_0__1);
		iv1__1_1__1_0_0_1 = IVector(7, aiv1__1_1__1_0_0_1);
		iv1__1_1__1_0_0_0 = IVector(7, aiv1__1_1__1_0_0_0);
		iv1_0_1_0_0_0_0 = IVector(7, aiv1_0_1_0_0_0_0);
		iv0_1_0_1_0_0_0 = IVector(7, aiv0_1_0_1_0_0_0);

		iv0_1_0_0_0_0__1 = IVector(7, aiv0_1_0_0_0_0__1);
		iv0_0_0_1_0_0_1 = IVector(7, aiv0_0_0_1_0_0_1);

		iv0_1_1_0_0_0_0 = IVector(7, aiv0_1_1_0_0_0_0);
		iv1_0_0_0_0_0__1 = IVector(7, aiv1_0_0_0_0_0__1);

		iv1_0_0_1_0_0_0 = IVector(7, aiv1_0_0_1_0_0_0);
		iv0_0_1_0_0_0_1 = IVector(7, aiv0_0_1_0_0_0_1);
	}

	IAffineTransformation iat_i_j, iat_im1_j, iat_i_jm1,
		iat_i_j_1, iat_i_1_j, iat_1_i_j,
		iat_im5_j, iat_im3_jm4,
		iat_i, iat_im5, iat_im6, iat_im1;
	void test_RDMGraph_initializeAffineTransformations() {
		iat_i_j = IAffineTransformation(4);
		iat_i_j.addTransformationRow(new IVector(iv1_0_0_0_0));
		iat_i_j.addTransformationRow(new IVector(iv0_1_0_0_0));

		iat_im1_j = IAffineTransformation(4);
		iat_im1_j.addTransformationRow(new IVector(iv1_0_0_0__1));
		iat_im1_j.addTransformationRow(new IVector(iv0_1_0_0_0));

		iat_i_jm1 = IAffineTransformation(4);
		iat_i_jm1.addTransformationRow(new IVector(iv1_0_0_0_0));
		iat_i_jm1.addTransformationRow(new IVector(iv0_1_0_0__1));

		iat_i_j_1 = IAffineTransformation(4);
		iat_i_j_1.addTransformationRow(new IVector(iv1_0_0_0_0));
		iat_i_j_1.addTransformationRow(new IVector(iv0_1_0_0_0));
		iat_i_j_1.addTransformationRow(new IVector(iv0_0_0_0_1));

		iat_i_1_j = IAffineTransformation(4);
		iat_i_1_j.addTransformationRow(new IVector(iv1_0_0_0_0));
		iat_i_1_j.addTransformationRow(new IVector(iv0_0_0_0_1));
		iat_i_1_j.addTransformationRow(new IVector(iv0_1_0_0_0));

		iat_1_i_j = IAffineTransformation(4);
		iat_1_i_j.addTransformationRow(new IVector(iv0_0_0_0_1));
		iat_1_i_j.addTransformationRow(new IVector(iv1_0_0_0_0));
		iat_1_i_j.addTransformationRow(new IVector(iv0_1_0_0_0));

		iat_im5_j = IAffineTransformation(4);
		iat_im5_j.addTransformationRow(new IVector(iv1_0_0_0__5));
		iat_im5_j.addTransformationRow(new IVector(iv0_1_0_0_0));

		iat_im3_jm4 = IAffineTransformation(4);
		iat_im3_jm4.addTransformationRow(new IVector(iv1_0_0_0__3));
		iat_im3_jm4.addTransformationRow(new IVector(iv0_1_0_0__4));

		iat_i = IAffineTransformation(3);
		iat_i.addTransformationRow(new IVector(iv1_0_0_0));

		iat_im5 = IAffineTransformation(3);
		iat_im5.addTransformationRow(new IVector(iv1_0_0__5));

		iat_im6 = IAffineTransformation(3);
		iat_im6.addTransformationRow(new IVector(iv1_0_0__6));

		iat_im1 = IAffineTransformation(3);
		iat_im1.addTransformationRow(new IVector(iv1_0_0__1));
	}

	IModuleRelation imr_imj_imj_m1, imr_imj_imj_1, imr_imj_imj,
						imr_i_i_0__j_j_0, imr_i_i_0__j_0_m1___0_j_1, imr_j_i_0__i_0_m1___0_j_1, imr_i_j_0__j_0_m1___0_i_1, imr_j_j_0__i_0_m1___0_i_1;
	void test_RDMGraph_initializeModuleRelations() {
		imr_imj_imj_m1 = IModuleRelation(2, 2, 3);
		imr_imj_imj_m1.addGeneratingVertex(new IVector(iv1__1_1__1_0_0__1));

		imr_imj_imj_1 = IModuleRelation(2, 2, 3);
		imr_imj_imj_1.addGeneratingVertex(new IVector(iv1__1_1__1_0_0_1));

		imr_imj_imj = IModuleRelation(2, 2, 3);
		imr_imj_imj.addGeneratingVertex(new IVector(iv1__1_1__1_0_0_0));

		imr_i_i_0__j_j_0 = IModuleRelation(2, 2, 3);
		imr_i_i_0__j_j_0.addGeneratingVertex(new IVector(iv1_0_1_0_0_0_0));
		imr_i_i_0__j_j_0.addGeneratingVertex(new IVector(iv0_1_0_1_0_0_0));

		imr_i_i_0__j_0_m1___0_j_1 = IModuleRelation(2, 2, 3);
		imr_i_i_0__j_0_m1___0_j_1.addGeneratingVertex(new IVector(iv1_0_1_0_0_0_0));
		imr_i_i_0__j_0_m1___0_j_1.addGeneratingVertex(new IVector(iv0_1_0_0_0_0__1));
		imr_i_i_0__j_0_m1___0_j_1.addGeneratingVertex(new IVector(iv0_0_0_1_0_0_1));

		imr_j_i_0__i_0_m1___0_j_1 = IModuleRelation(2, 2, 3);
		imr_j_i_0__i_0_m1___0_j_1.addGeneratingVertex(new IVector(iv0_1_1_0_0_0_0));
		imr_j_i_0__i_0_m1___0_j_1.addGeneratingVertex(new IVector(iv1_0_0_0_0_0__1));
		imr_j_i_0__i_0_m1___0_j_1.addGeneratingVertex(new IVector(iv0_0_0_1_0_0_1));

		imr_i_j_0__j_0_m1___0_i_1 = IModuleRelation(2, 2, 3);
		imr_i_j_0__j_0_m1___0_i_1.addGeneratingVertex(new IVector(iv1_0_0_1_0_0_0));
		imr_i_j_0__j_0_m1___0_i_1.addGeneratingVertex(new IVector(iv0_1_0_0_0_0__1));
		imr_i_j_0__j_0_m1___0_i_1.addGeneratingVertex(new IVector(iv0_0_1_0_0_0_1));

		imr_j_j_0__i_0_m1___0_i_1 = IModuleRelation(2, 2, 3);
		imr_j_j_0__i_0_m1___0_i_1.addGeneratingVertex(new IVector(iv0_1_0_1_0_0_0));
		imr_j_j_0__i_0_m1___0_i_1.addGeneratingVertex(new IVector(iv1_0_0_0_0_0__1));
		imr_j_j_0__i_0_m1___0_i_1.addGeneratingVertex(new IVector(iv0_0_1_0_0_0_1));
	}

	typedef CVertex<CPolyheder<I>, IAffineTransformation> RDMVertex;
	typedef CEdge<CPolyheder<I>, IAffineTransformation> RDMEdge;
	typedef CRDMGraph<I> RDMGraph;

	typedef CVertex<CPolyheder<I>, IModuleRelation> RADVertex;
	typedef CEdge<CPolyheder<I>, IModuleRelation> RADEdge;
	typedef CRADGraph<I, IModuleRelation> RADGraph;

	typedef CVertex<CPolyheder<I>, ILatticeRelation> RADLatticeVertex;
	typedef CEdge<CPolyheder<I>, ILatticeRelation> RADLatticeEdge;
	typedef CRADGraph<I, ILatticeRelation> RADLatticeGraph;

	typedef CVertex<CPolyheder<I>, IModuleRelation> RADModuleVertex;
	typedef CEdge<CPolyheder<I>, IModuleRelation> RADModuleEdge;
	typedef CRADGraph<I, IModuleRelation> RADModuleGraph;

	RDMGraph g;
	vector<RDMVertex* > cV;
	vector<RDMVertex* > dV;

	void test_RDMGraph_SetupGraph(int compVertices, int dataVertices, int paramCount) {
		test_RDMGraph_initializeIVectors();
		test_RDMGraph_initializeAffineTransformations();
		test_RDMGraph_initializeModuleRelations();

		g = RDMGraph();
		cV = vector<RDMVertex* >(compVertices);
		for (unsigned int q = 0; q < cV.size(); q++) cV[q] = &g.addComputationVertex(IPolyheder());

		dV = vector<RDMVertex* >(dataVertices);
		for (unsigned int q = 0; q < dV.size(); q++) dV[q] = &g.addDataVertex(IPolyheder());

		g.parmCount = paramCount;

		ASSERT(g.getVertexCount() == compVertices + dataVertices);
	}

	RADGraph* radg;
	RADLatticeGraph* radgLattice;

	vector<RADVertex* > w;
	vector<RADLatticeVertex* > wLattice;

	vector<vector<vector<RADEdge*> > > resultEdges;
	vector<vector<RADLatticeEdge*> > resultEdgesLattice;

	vector<ILattice*> distanceLattices;

	vector<vector<vector<RADEdge*> > > extractResultEdges(RADGraph* gFW, vector<RADVertex* >& w) {
		vector<vector<vector<RADEdge*> > > resultEdges
			= vector<vector<vector<RADEdge*> > >(w.size(), vector<vector<RADEdge*> >(w.size()));

		for (unsigned int q = 0; q < w.size(); q++) {
			for (unsigned int r = 0; r < w.size(); r++) {
				resultEdges[q][r] = gFW->getEdgesBetween(*w[q], *w[r]);
			}
		}

		return resultEdges;
	}

	vector<vector<RADLatticeEdge*> > extractResultLatticeEdges(RADLatticeGraph* gFW, vector<RADLatticeVertex* >& w) {
		vector<vector<RADLatticeEdge*> > resultEdgesLattice
			= vector<vector<RADLatticeEdge*> >(w.size(), vector<RADLatticeEdge*>(w.size()));

		for (unsigned int q = 0; q < w.size(); q++) {
			for (unsigned int r = 0; r < w.size(); r++) {
				vector<RADLatticeEdge*> edges = gFW->getEdgesBetween(*w[q], *w[r]);
				ASSERT((edges.size() == 0) || (edges.size() == 1));
				resultEdgesLattice[q][r] = (edges.size() == 1) ? edges[0] : NULL;
			}
		}

		return resultEdgesLattice;
	}

	void determineDistanceLattices(ofstream* ofsHtml = NULL) {
		distanceLattices = vector<ILattice*>(resultEdgesLattice.size());
		if (ofsHtml != NULL ) *ofsHtml << "<br/><br/><b>Determining distance lattices:</b><br/>";
		for (unsigned int q = 0; q < resultEdgesLattice.size(); q++) {
			if (resultEdgesLattice[q][q] != NULL) {
				if (ofsHtml != NULL ) *ofsHtml << "Statement " << q << ":<br/>";
				ILatticeRelation& lr = resultEdgesLattice[q][q]->getData();

				distanceLattices[q] = new ILattice(calculateDistanceLattice(lr));
				//*ofsHtml << "Generating matrix of projected result:<br/>";
				if (ofsHtml != NULL ) *ofsHtml << distanceLattices[q]->toStringHtml();
				if (ofsHtml != NULL ) *ofsHtml << "<br/>";
			} else distanceLattices[q] = NULL;
		}
	}

	void performSpacePartitioning(ostream* osLatex = NULL, ostream* osHtml = NULL) {
		if (osHtml != NULL) *osHtml << "<b>Starting Affine Space Partitioning</b><br/>";

		CompGraphSpacePartitioner<I> sp = CompGraphSpacePartitioner<I>(&g);
		sp.osHtml = osHtml;
		sp.performAffineSetSpacePartitioning();
		radg = sp.dependenceFlatClosure;

		w = vector<RADVertex* >();
		for (int q = 0; q < g.getComputationVertexCount(); q++) w.push_back(&radg->getVertex(q));

		resultEdges = extractResultEdges(radg, w);
	}

	struct FWResultGraphData {
		vector<string> statementNames;
	};

	void encodeFWResultsInHtml(ostream& ofs, CVector<int> lowerBounds, CVector<int> upperBounds,
			                   FWResultGraphData* fw = NULL) {
		int basicWidth = 50000;
		CSVGEElementCollection table = CSVGEElementCollection();
		int spaceW = basicWidth/(cV.size()*(spaceFact + 1) - 1);
		int cellW = spaceFact*spaceW;
		for (unsigned int q = 0; q < cV.size(); q++) {
			for (unsigned int r = 0; r < cV.size(); r++) {
				if (resultEdgesLattice[q][r] != NULL) {
					ILattice projLattice = resultEdgesLattice[q][r]->getData();

					//ofs << resultEdgesLattice[q][r]->getData().toStringHtml();
					CSVGELattice<I>* lat = new CSVGELattice<I>(projLattice, lowerBounds, upperBounds, basicWidth, true);
					lat->horAxisTitle = "Statement " + intToStr(r);
					lat->verAxisTitle = "Statement " + intToStr(q);
					lat->swapped = true;

					table.elements.push_back(new CSVGEViewBox(lat,
							r*(cellW + spaceW), q*(cellW + spaceW),
							cellW, cellW, basicWidth, basicWidth));
				}
			}
		}

		CSVG(&table, 700, 700).encode(ofs);
	}

	void encodeFWResultsInPdf(string pdfFile, CVector<int> lowerBounds, CVector<int> upperBounds,
			                   FWResultGraphData* fw = NULL) {
		int basicWidth = 50000;
		CSVGEElementCollection table = CSVGEElementCollection();
		int spaceW = basicWidth/(cV.size()*(spaceFact + 1) - 1);
		int cellW = spaceFact*spaceW;
		for (unsigned int q = 0; q < cV.size(); q++) {
			for (unsigned int r = 0; r < cV.size(); r++) {
				if (resultEdgesLattice[q][r] != NULL) {
					ILattice projLattice = resultEdgesLattice[q][r]->getData();

					//ofs << resultEdgesLattice[q][r]->getData().toStringHtml();
					CSVGELattice<I>* lat = new CSVGELattice<I>(projLattice, lowerBounds, upperBounds, basicWidth, true);
					string ij = "ij";
					lat->horAxisTitle = ij[r];
					lat->verAxisTitle = ij[q];
					lat->swapped = true;

					table.elements.push_back(new CSVGEViewBox(lat,
							r*(cellW + spaceW), q*(cellW + spaceW),
							cellW, cellW, basicWidth, basicWidth));
				}
			}
		}

		CSVG(&table, 700, 700).saveToPdf(pdfFile);
	}

	CSVGElement* elementsToTable(vector<CSVGElement*> elements, int itemsPerRow) {
		int basicWidth = 50000;
		CSVGEElementCollection* table = new CSVGEElementCollection();
		int spaceW = basicWidth/(itemsPerRow*(spaceFact + 1) - 1);
		int cellW = spaceFact*spaceW;
		for (unsigned int elIx = 0; elIx < elements.size(); ) {
			for (int q = 0; (q < itemsPerRow) && (elIx < elements.size()); q++) {
				table->elements.push_back(new CSVGEViewBox(elements[elIx],
							q*(cellW + spaceW), (elIx / itemsPerRow)*(cellW + spaceW),
							cellW, cellW, basicWidth, basicWidth));

				elIx++;
			}
		}

		return table;
	}

	vector<CSVGElement*> elementsToTables(vector<CSVGElement*> elements, int itemsPerRow, int itemsPerTable) {
		vector<CSVGElement*> result;

		for (unsigned int elIx = 0; elIx < elements.size(); ) {
			vector<CSVGElement*> tableElements;
			for (int q = 0; (q < itemsPerTable) && (elIx < elements.size()); q++) {
				tableElements.push_back(elements[elIx++]);
			}
			result.push_back(elementsToTable(tableElements, itemsPerRow));
		}

		return result;
	}

	string verboseOutputPrefix, verboseOutputPath, fwInputPdf, initialLatticesPdf;
	CVector<int> lowerBounds, upperBounds;

	void performLatticeSpacePartitioning(ostream* osLatex = NULL, ostream* osHtml = NULL) {
		g.osLatex = osLatex;
		g.osHtml = osHtml;
		g.verboseOutputPath = verboseOutputPath;
		g.verboseOutputPrefix = verboseOutputPrefix;
		g.latticeLowerBounds = lowerBounds;
		g.latticeUpperBounds = upperBounds;
		g.latticeBasicWidth = 50000;

		if (osHtml != NULL) *osHtml << "<b>Starting Lattice-Based Space Partitioning</b><br/>";

		CompGraphSpacePartitioner<I> sp = CompGraphSpacePartitioner<I>(&g);
		sp.osHtml = osHtml;
		sp.performAffineLatticeSpacePartitioning();

		/*g.osHtml = NULL;
		RADLatticeGraph flApproxRSDG = g.calculateFlatLatticeApproximatedRSDG();

		if (osLatex != NULL) {
			vector<CSVGElement* > initialEdges = g.resultingSVGLattices;
			vector<CSVGElement*> tables = elementsToTables(initialEdges, 2, 2);
			for (unsigned int q = 0; q < tables.size(); q++) {
				CSVG(tables[q], 700, 700).saveToPdf(initialLatticesPdf + intToStr(q) + ".pdf");
				delete tables[q];
			}

			RADLatticeGraph* combinedFLARG = flApproxRSDG.calcApproximateDisjunctiveCombinationReduction();
			wLattice = vector<RADLatticeVertex* >();
			for (unsigned int q = 0; q < cV.size(); q++) wLattice.push_back(&combinedFLARG->getVertex(q));
			resultEdgesLattice = extractResultLatticeEdges(combinedFLARG, wLattice);

			encodeFWResultsInPdf(fwInputPdf, lowerBounds, upperBounds);

			delete combinedFLARG;
		}*/

		radgLattice = sp.dependenceLatticeClosure;

		wLattice = vector<RADLatticeVertex* >();
		for (unsigned int q = 0; q < cV.size(); q++) wLattice.push_back(&radgLattice->getVertex(q));
		resultEdgesLattice = extractResultLatticeEdges(radgLattice, wLattice);
	}

	/* Trivial 2D parallelism
	 * for ( i = 0 ; i <= n ; i++) {
	 * 	 for ( j = 0 ; j <= m ; j++) {
	 *     A[i][j] := A[i][j] + B[i][j];
	 *     B[i][j] := A[i][j] * B[i][j];
	 *   }
	 * }
	 */
	void CRDMGraph_SpacePartition_Trivial2DParallelism() {
		test_RDMGraph_SetupGraph(2, 2, 2);

		cV[0]->setData(IPolyheder(4));
		cV[1]->setData(IPolyheder(4));
		dV[0]->setData(IPolyheder(2));
		dV[1]->setData(IPolyheder(2));

		g.addEdge(cV[0], dV[0], iat_i_j);
		g.addEdge(cV[0], dV[1], iat_i_j);
		g.addEdge(cV[1], dV[0], iat_i_j);
		g.addEdge(cV[1], dV[1], iat_i_j);

		performSpacePartitioning();

		ASSERT((resultEdges[0][1].size() == 1)
				&& (resultEdges[0][1][0]->getData() == IMR_Con(2, 2, 3, 2, IV_(7,  1,  0, -1,  0,  0,  0,  0),
						                                              IV_(7,  0,  1,  0, -1,  0,  0,  0))));
		ASSERT((resultEdges[1][0].size() == 1)
				&& (resultEdges[1][0][0]->getData() == IMR_Con(2, 2, 3, 2, IV_(7,  1,  0, -1,  0,  0,  0,  0),
																	  IV_(7,  0,  1,  0, -1,  0,  0,  0))));
		ASSERT((resultEdges[0][0].size() == 1)
				&& (resultEdges[0][0][0]->getData() == IMR_Con(2, 2, 3, 2, IV_(7,  1,  0, -1,  0,  0,  0,  0),
						   											  IV_(7,  0,  1,  0, -1,  0,  0,  0))));
		ASSERT((resultEdges[1][1].size() == 1)
				&& (resultEdges[1][1][0]->getData() == IMR_Con(2, 2, 3, 2, IV_(7,  1,  0, -1,  0,  0,  0,  0),
							  										  IV_(7,  0,  1,  0, -1,  0,  0,  0))));
	}

	/* Default Lim & Lam example
	 * for ( i = 0 ; i <= n ; i++) {
	 * 	 for ( j = 0 ; j <= m ; j++) {
	 *     A[i][j] := A[i][j] + B[i - 1][j];
	 *     B[i][j] := A[i][j - 1] * B[i][j];
	 *   }
	 * }
	 */
	void CRDMGraph_AffineSet_SpacePartition_LimLam() {
		test_RDMGraph_SetupGraph(2, 2, 2);

		cV[0]->setData(IPolyheder(4));
		cV[1]->setData(IPolyheder(4));
		dV[0]->setData(IPolyheder(2));
		dV[1]->setData(IPolyheder(2));

		g.addEdge(cV[0], dV[0], iat_i_j);
		g.addEdge(cV[0], dV[1], iat_im1_j);
		g.addEdge(cV[1], dV[0], iat_i_jm1);
		g.addEdge(cV[1], dV[1], iat_i_j);

		performSpacePartitioning();

		ASSERT((resultEdges[0][1].size() == 1)
				&& (resultEdges[0][1][0]->getData() == IMR_Con(2, 2, 3, 1, IV_(7,  1, -1, -1,   1,  0,  0, -1)))); // imr_imj_imj_m1
		ASSERT((resultEdges[1][0].size() == 1)
				&& (resultEdges[1][0][0]->getData() == IMR_Con(2, 2, 3, 1, IV_(7,  1, -1, -1,   1,  0,  0,  1)))); // imr_imj_imj_1
		ASSERT((resultEdges[0][0].size() == 1)
				&& (resultEdges[0][0][0]->getData() == IMR_Con(2, 2, 3, 1, IV_(7,  1, -1, -1,   1,  0,  0,  0)))); // imr_imj_imj
		ASSERT((resultEdges[1][1].size() == 1)
				&& (resultEdges[1][1][0]->getData() == IMR_Con(2, 2, 3, 1, IV_(7,  1, -1, -1,   1,  0,  0,  0)))); // imr_imj_imj
	}

	/* Default Lim & Lam example
	 * for ( i = 0 ; i <= n ; i++) {
	 * 	 for ( j = 0 ; j <= m ; j++) {
	 *     A[i][j] := A[i][j] + B[i - 1][j];
	 *     B[i][j] := A[i][j - 1] * B[i][j];
	 *   }
	 * }
	 */
	void CRDMGraph_AffineLattice_SpacePartition_LimLam() {
		test_RDMGraph_SetupGraph(2, 2, 2);

		cV[0]->setData(IPolyheder(4));
		cV[1]->setData(IPolyheder(4));
		dV[0]->setData(IPolyheder(2));
		dV[1]->setData(IPolyheder(2));

		g.addEdge(cV[0], dV[0], iat_i_j);
		g.addEdge(cV[0], dV[1], iat_im1_j);
		g.addEdge(cV[1], dV[0], iat_i_jm1);
		g.addEdge(cV[1], dV[1], iat_i_j);

		performLatticeSpacePartitioning();
	}

	void encodeDistanceLatticesInHtml(ostream& ofs, CVector<int> lowerBounds, CVector<int> upperBounds) {
		int basicWidth = 50000;
		CSVGEElementCollection table = CSVGEElementCollection();
		int spaceW = basicWidth/(cV.size()*(spaceFact + 1) - 1);
		int cellW = spaceFact*spaceW;
		for (unsigned int q = 0; q < cV.size(); q++) {
			ILattice* dL = distanceLattices[q];

			//ofs << dL->toStringHtml();
			CSVGELattice<I>* lat = new CSVGELattice<I>(*dL, lowerBounds, upperBounds, basicWidth, false);
			table.elements.push_back(new CSVGEViewBox(lat,
					q*(cellW + spaceW), 0,
					cellW, cellW, basicWidth, basicWidth));
			string ij = "ij";
			lat->horAxisTitle = ij[q];
		}

		CSVG(&table, 300, 300).encode(ofs);
	}

	void encodeDistanceLatticesInPdf(string pdfName, CVector<int> lowerBounds, CVector<int> upperBounds) {
		int basicWidth = 50000;
		CSVGEElementCollection table = CSVGEElementCollection();
		int spaceW = basicWidth/(cV.size()*(spaceFact + 1) - 1);
		int cellW = spaceFact*spaceW;
		for (unsigned int q = 0; q < cV.size(); q++) {
			ILattice* dL = distanceLattices[q];

			//ofs << dL->toStringHtml();
			CSVGELattice<I>* lat = new CSVGELattice<I>(*dL, lowerBounds, upperBounds, basicWidth, false);
			table.elements.push_back(new CSVGEViewBox(lat,
					q*(cellW + spaceW), 0,
					cellW, cellW, basicWidth, basicWidth));
			string ij = "ij";
			lat->horAxisTitle = ij[q];
		}

		CSVG(&table, 300, 300).saveToPdf(pdfName);
	}

	/* Simple 1D non-uniform loopnest for lattice-based parallelism
	 * for ( i = 0 ; i <= 100 ; i++) {
	 *   A[i] := A[i] + A[i - 6];
	 * }
	 */
	void CRDMGraph_SpacePartition_1DUniformLattice_SingleStatement() {
		test_RDMGraph_SetupGraph(1, 1, 0);

		cV[0]->setData(IPolyheder(1));
		dV[0]->setData(IPolyheder(1));

		g.addEdge(cV[0], dV[0], IAT(1, IV_(2,  1,  0)));
		g.addEdge(cV[0], dV[0], IAT(1, IV_(2,  1, -6)));

		ofstream ofs("examples/CRDMGraph_SpacePartition_1DUniformLattice_SingleStatement.xhtml");
		ofs << "<html xmlns=\"http://www.w3.org/1999/xhtml\"><body>";

		performLatticeSpacePartitioning(NULL, &ofs);

		ASSERT(resultEdgesLattice[0][0]->getData() == ILR(1, 1, 1, 3, IV_(3, 1, 1, 0),
																	  IV_(3, 0, 6, 0),
																	  IV_(3, 0, 0, 1)));

		encodeDistanceLatticesInHtml(ofs, IntV(1, -7), IntV(1, 7));

		ofs << "</body></html>"; ofs.close();
	}

	/* Simple 1D non-uniform loopnest for lattice-based parallelism,
	 * but it should yield flat lattice parallelism
	 * for ( i = 0 ; i <= 100 ; i++) {
	 *   A[i] := A[i - 3] + A[4*i];
	 * }
	 */
	void CRDMGraph_SpacePartition_1DUniformLattice_SingleStatement_NonUniform() {
		test_RDMGraph_SetupGraph(1, 1, 0);

		cV[0]->setData(IP(2, IV_(2,  1,  0),
							 IV_(2, -1, 100)));
		dV[0]->setData(IPolyheder(1));

		g.addEdge(cV[0], dV[0], IAT(1, IV_(2,  4,  0)));
		g.addEdge(cV[0], dV[0], IAT(1, IV_(2,  1, -3)));
		g.addEdge(cV[0], dV[0], IAT(1, IV_(2,  1,  0)));

		performLatticeSpacePartitioning();

		ofstream ofs("CRDMGraph_SpacePartition_1DUniformLattice_SingleStatement_NonUniform.xhtml");
		ofs << "<html xmlns=\"http://www.w3.org/1999/xhtml\"><body>";
		encodeFWResultsInHtml(ofs, IntV(2, 0, 0), IntV(2, 15, 15));

		ASSERT(resultEdgesLattice[0][0]->getData() == ILR(1, 1, 1, 3, IV_(3, 1, 1, 0),
																	  IV_(3, 0, 3, 0),
																	  IV_(3, 0, 0, 1)));

		ofs << "<hr />";
		determineDistanceLattices();

		encodeDistanceLatticesInHtml(ofs, IntV(1, -7), IntV(1, 7));

		ofs << "</body></html>"; ofs.close();
	}

	void CRDMGraph_AffineLattice_SpacePartition_1DUniformLattice_NonUniform() {
		test_RDMGraph_SetupGraph(2, 2, 0);

		cV[0]->setData(IP(2, IV_(2,  1,  0),
							 IV_(2, -1, 10)));
		dV[0]->setData(IPolyheder(1));
		cV[1]->setData(IP(2, IV_(2,  1,  0),
							 IV_(2, -1, 10)));
		dV[1]->setData(IPolyheder(1));

		g.addEdge(cV[0], dV[0], IAT(1, IV_(2,  1,  1)));
		g.addEdge(cV[0], dV[0], IAT(1, IV_(2,  1, -2)));
		g.addEdge(cV[1], dV[0], IAT(1, IV_(2,  4,  0)));

		g.addEdge(cV[0], dV[1], IAT(1, IV_(2,  5,  0)));
		g.addEdge(cV[1], dV[1], IAT(1, IV_(2,  5,  4)));

		ofstream ofs("examples/CRDMGraph_AffineLattice_SpacePartition_1DUniformLattice_NonUniform.xhtml");
		ofs << "<html xmlns=\"http://www.w3.org/1999/xhtml\"><body>";
		ofs << "<font size=\"+3\"><pre>"
		    << "for (i := 0 to 10) {\n"
		    << "  A[i + 1] = A[i - 2] * B[5*i];\n"
		    << "  B[5*i + 4] = A[4*i] + B[5*i + 4];\n"
		    << "}\n"
		    << "</pre></font>";

		performLatticeSpacePartitioning(NULL, &ofs);

		/*ASSERT(resultEdgesLattice[0][0]->getData() == ILR(1, 1, 1, 3, IV_(3, 1, 1, 0),
																 	  IV_(3, 0, 3, 0),
																	  IV_(3, 0, 0, 1)));
		ASSERT(resultEdgesLattice[0][1]->getData() == ILR(1, 1, 1, 3, IV_(3, 1, 0, 2),
						  											  IV_(3, 0, 1, 1),
						  											  IV_(3, 0, 0, 3)));
		ASSERT(resultEdgesLattice[0][1]->getData() == ILR(1, 1, 1, 3, IV_(3, 1, 0, 1),
						  											  IV_(3, 0, 1, 2),
						  											  IV_(3, 0, 0, 3)));
		ASSERT(resultEdgesLattice[1][1]->getData() == ILR(1, 1, 1, 3, IV_(3, 1, 1, 0),
																	  IV_(3, 0, 3, 0),
																	  IV_(3, 0, 0, 1)));*/

		ofs << "</body></html>"; ofs.close(); //ofsL.close();
	}

	void CRDMGraph_AffineSet_SpacePartition_1DUniformLattice_NonUniform() {
		test_RDMGraph_SetupGraph(2, 2, 0);

		cV[0]->setData(IP(2, IV_(2,  1,  0),
							 IV_(2, -1, 10)));
		dV[0]->setData(IPolyheder(1));
		cV[1]->setData(IP(2, IV_(2,  1,  0),
							 IV_(2, -1, 10)));
		dV[1]->setData(IPolyheder(1));

		g.addEdge(cV[0], dV[0], IAT(1, IV_(2,  1,  1)));
		g.addEdge(cV[0], dV[0], IAT(1, IV_(2,  1, -2)));
		g.addEdge(cV[1], dV[0], IAT(1, IV_(2,  4,  0)));

		g.addEdge(cV[0], dV[1], IAT(1, IV_(2,  5,  0)));
		g.addEdge(cV[1], dV[1], IAT(1, IV_(2,  5,  4)));

		ofstream ofs("examples/CRDMGraph_AffineSet_SpacePartition_1DUniformLattice_NonUniform.xhtml");
		ofs << "<html xmlns=\"http://www.w3.org/1999/xhtml\"><body>";
		ofs << "<font size=\"+3\"><pre>"
		    << "for (i := 0 to 10) {\n"
		    << "  A[i + 1] = A[i - 2] * B[5*i];\n"
		    << "  B[5*i + 4] = A[4*i] + B[5*i + 4];\n"
		    << "}\n"
		    << "</pre></font>";

		performSpacePartitioning(NULL, &ofs);

		/*ASSERT(resultEdgesLattice[0][0]->getData() == ILR(1, 1, 1, 3, IV_(3, 1, 1, 0),
																 	  IV_(3, 0, 3, 0),
																	  IV_(3, 0, 0, 1)));
		ASSERT(resultEdgesLattice[0][1]->getData() == ILR(1, 1, 1, 3, IV_(3, 1, 0, 2),
						  											  IV_(3, 0, 1, 1),
						  											  IV_(3, 0, 0, 3)));
		ASSERT(resultEdgesLattice[0][1]->getData() == ILR(1, 1, 1, 3, IV_(3, 1, 0, 1),
						  											  IV_(3, 0, 1, 2),
						  											  IV_(3, 0, 0, 3)));
		ASSERT(resultEdgesLattice[1][1]->getData() == ILR(1, 1, 1, 3, IV_(3, 1, 1, 0),
																	  IV_(3, 0, 3, 0),
																	  IV_(3, 0, 0, 1)));*/

		ofs << "</body></html>"; ofs.close(); //ofsL.close();
	}

	/*
	 for (i = 0; i <= 10; i++) {
	   A[2*i - 4] := A[2*i + 8] + B[i]
	   B[i - 1] := B[i + 3] + A[2*i - 6]
	 }
	*/
	void CRDMGraph_SpacePartition_1DUniformLattice_NonUniform2() {
		test_RDMGraph_SetupGraph(2, 2, 0);

		cV[0]->setData(IPolyheder(1));
		dV[0]->setData(IPolyheder(1));
		cV[1]->setData(IPolyheder(1));
		dV[1]->setData(IPolyheder(1));

		g.addEdge(cV[0], dV[0], IAT(1, IV_(2,  2, -4)));
		g.addEdge(cV[0], dV[0], IAT(1, IV_(2,  2,  8)));
		g.addEdge(cV[0], dV[1], IAT(1, IV_(2,  1,  0)));

		g.addEdge(cV[1], dV[0], IAT(1, IV_(2,  2, -6)));
		g.addEdge(cV[1], dV[1], IAT(1, IV_(2,  1, -1)));
		g.addEdge(cV[1], dV[1], IAT(1, IV_(2,  1,  8)));

		ofstream ofsL("../CommunicationFreeParallelisation_Conference/src/Examples/example1.tex");

		ofstream ofs("examples/CRDMGraph_SpacePartition_1DUniformLattice_NonUniform2.xhtml");
		//ofs << "<html xmlns=\"http://www.w3.org/1999/xhtml\"><body>";
		performLatticeSpacePartitioning(&ofsL, &ofs);


		encodeFWResultsInHtml(ofs, IntV(2, 0, 0), IntV(2, 10, 10));

		/*ASSERT(resultEdgesLattice[0][0]->getData() == ILR(1, 1, 1, 3, IV_(3, 1, 1, 0),
																 	  IV_(3, 0, 3, 0),
																	  IV_(3, 0, 0, 1)));
		ASSERT(resultEdgesLattice[0][1]->getData() == ILR(1, 1, 1, 3, IV_(3, 1, 0, 2),
						  											  IV_(3, 0, 1, 1),
						  											  IV_(3, 0, 0, 3)));
		ASSERT(resultEdgesLattice[0][1]->getData() == ILR(1, 1, 1, 3, IV_(3, 1, 0, 1),
						  											  IV_(3, 0, 1, 2),
						  											  IV_(3, 0, 0, 3)));
		ASSERT(resultEdgesLattice[1][1]->getData() == ILR(1, 1, 1, 3, IV_(3, 1, 1, 0),
																	  IV_(3, 0, 3, 0),
																	  IV_(3, 0, 0, 1)));*/

		ofs << "<hr />";
		determineDistanceLattices();

		encodeDistanceLatticesInHtml(ofs, IntV(1, -7), IntV(1, 7));

		ofs << "</body></html>"; ofs.close(); ofsL.close();
	}

	/*
	 for (i = 0; i <= 10; i++) {
	   A[2*i - 4] := A[2*i + 8] + B[i]
	   B[i - 1] := B[i + 8] + A[2*i - 6]
	 }
	*/
	void CRDMGraph_SpacePartition_1DUniformLattice_NonUniform3() {
		test_RDMGraph_SetupGraph(2, 2, 0);

		cV[0]->setData(IP(2, IV_(2,  1,  0),
				             IV_(2, -1, 10)));
		dV[0]->setData(IPolyheder(1));
		cV[1]->setData(IP(2, IV_(2,  1,  0),
				             IV_(2, -1, 10)));
		dV[1]->setData(IPolyheder(1));

		g.addEdge(cV[0], dV[0], IAT(1, IV_(2,  2, -4)));
		g.addEdge(cV[0], dV[0], IAT(1, IV_(2,  2,  8)));
		g.addEdge(cV[0], dV[1], IAT(1, IV_(2,  1,  0)));

		g.addEdge(cV[1], dV[0], IAT(1, IV_(2,  2, -6)));
		g.addEdge(cV[1], dV[1], IAT(1, IV_(2,  1, -1)));
		g.addEdge(cV[1], dV[1], IAT(1, IV_(2,  1,  8)));

		ofstream ofsL("../CommunicationFreeParallelism_Conference/src/Examples/example1.tex");

		verboseOutputPath = "../CommunicationFreeParallelism_Conference/src/";
		verboseOutputPrefix = "Examples/x1/";
		lowerBounds = IntV(2, 0, 0);
		upperBounds = IntV(2, 10, 10);
		fwInputPdf = "../CommunicationFreeParallelism_Conference/src/Examples/x1/FWInput.pdf";
		initialLatticesPdf = "../CommunicationFreeParallelism_Conference/src/Examples/x1/initialLattices";


		ofstream ofs("examples/CRDMGraph_SpacePartition_1DUniformLattice_NonUniform3.xhtml");
		ofs << "<html xmlns=\"http://www.w3.org/1999/xhtml\"><body>";
		performLatticeSpacePartitioning(&ofsL, &ofs);

		//encodeFWResultsInHtml(ofs, IntV(2, 0, 0), IntV(2, 10, 10));

		encodeFWResultsInPdf("../CommunicationFreeParallelism_Conference/src/Examples/x1/FWResult.pdf", IntV(2, 0, 0), IntV(2, 10, 10));
		/*ASSERT(resultEdgesLattice[0][0]->getData() == ILR(1, 1, 1, 3, IV_(3, 1, 1, 0),
																 	  IV_(3, 0, 3, 0),
																	  IV_(3, 0, 0, 1)));
		ASSERT(resultEdgesLattice[0][1]->getData() == ILR(1, 1, 1, 3, IV_(3, 1, 0, 2),
						  											  IV_(3, 0, 1, 1),
						  											  IV_(3, 0, 0, 3)));
		ASSERT(resultEdgesLattice[0][1]->getData() == ILR(1, 1, 1, 3, IV_(3, 1, 0, 1),
						  											  IV_(3, 0, 1, 2),
						  											  IV_(3, 0, 0, 3)));
		ASSERT(resultEdgesLattice[1][1]->getData() == ILR(1, 1, 1, 3, IV_(3, 1, 1, 0),
																	  IV_(3, 0, 3, 0),
																	  IV_(3, 0, 0, 1)));*/

		ofs << "<hr />";
		determineDistanceLattices();

		//encodeDistanceLatticesInHtml(ofs, IntV(1, -7), IntV(1, 7));

		ofstream ofsDist("examples/CRDMGraph_SpacePartition_1DUniformLattice_NonUniform3Distance.svg");
		//encodeDistanceLatticesInHtml(ofsDist, IntV(1, -7), IntV(1, 7));
		ofsDist.close();

		//encodeDistanceLatticesInPdf("../CommunicationFreeParallelism_Conference/src/Examples/x1/DistLat.pdf", IntV(1, -7), IntV(1, 7));

		ofs << "</body></html>"; ofs.close(); ofsL.close();
	}

	void CRDMGraph_SpacePartition_1DUniformLattice_NonUniform3_AffineSet() {
		test_RDMGraph_SetupGraph(2, 2, 0);

		cV[0]->setData(IP(2, IV_(2,  1,  0),
				             IV_(2, -1, 10)));
		dV[0]->setData(IPolyheder(1));
		cV[1]->setData(IP(2, IV_(2,  1,  0),
				             IV_(2, -1, 10)));
		dV[1]->setData(IPolyheder(1));

		g.addEdge(cV[0], dV[0], IAT(1, IV_(2,  2, -4)));
		g.addEdge(cV[0], dV[0], IAT(1, IV_(2,  2,  8)));
		g.addEdge(cV[0], dV[1], IAT(1, IV_(2,  1,  0)));

		g.addEdge(cV[1], dV[0], IAT(1, IV_(2,  2, -6)));
		g.addEdge(cV[1], dV[1], IAT(1, IV_(2,  1, -1)));
		g.addEdge(cV[1], dV[1], IAT(1, IV_(2,  1,  8)));

		lowerBounds = IntV(2, 0, 0);
		upperBounds = IntV(2, 10, 10);
		ofstream ofs("examples/CRDMGraph_SpacePartition_1DUniformLattice_NonUniform3_AffineSet.xhtml");
		ofs << "<html xmlns=\"http://www.w3.org/1999/xhtml\"><body>";
		performSpacePartitioning(NULL, &ofs);

		//encodeFWResultsInHtml(ofs, IntV(2, 0, 0), IntV(2, 10, 10));


		ofs << "</body></html>"; ofs.close();
	}

	/*
	 * for i := 0 to n do begin
	 * 	 for j := 0 to n do begin
	 *     A[i][j] := A[i + 2][j] + A[i][j + 2];
	 *   end;
	 * end;
	 */
	void CRDMGraph_SpacePartition_Param() {
		test_RDMGraph_SetupGraph(1, 1, 1);

		cV[0]->setData(IPolyheder(1));
		dV[0]->setData(IPolyheder(1));

		g.addEdge(cV[0], dV[0], IAT(1, IV_(3,  1,  0,  0)));
		g.addEdge(cV[0], dV[0], IAT(1, IV_(3,  1,  1,  0)));

		ofstream ofs("examples/CRDMGraph_SpacePartition_Param.xhtml");
		ofs << "<html xmlns=\"http://www.w3.org/1999/xhtml\"><body>";

		performLatticeSpacePartitioning(NULL, &ofs);

		//encodeDistanceLatticesInHtml(ofs, IntV(2, -7, -7), IntV(2, -7, -7));

		ofs << "</body></html>"; ofs.close();
	}

	/*
	 * for i := 0 to n do begin
	 * 	 for j := 0 to n do begin
	 *     A[i][j] := A[i + 2][j] + A[i][j + 2];
	 *   end;
	 * end;
	 */
	void CRDMGraph_SpacePartition_PACTReview_AffineLattice_Param() {
		test_RDMGraph_SetupGraph(1, 1, 1);

		cV[0]->setData(IP(4, IV_(4,  1,  0,  0,  0),
				             IV_(4, -1,  0,  1,  0),
				             IV_(4,  0,  1,  0,  0),
				             IV_(4,  0, -1,  1,  0)));
		dV[0]->setData(IPolyheder(2));

		g.addEdge(cV[0], dV[0], IAT(2, IV_(4,  1,  0,  0,  0),
				                       IV_(4,  0,  1,  0,  0)));
		g.addEdge(cV[0], dV[0], IAT(2, IV_(4,  1,  0,  0,  2),
				                       IV_(4,  0,  1,  0,  0)));
		g.addEdge(cV[0], dV[0], IAT(2, IV_(4,  1,  0,  0,  0),
				                       IV_(4,  0,  1,  0,  2)));

		ofstream ofs("examples/CRDMGraph_SpacePartition_PACTReview_AffineLattice_Param.xhtml");
		ofs << "<html xmlns=\"http://www.w3.org/1999/xhtml\"><body>";

		performLatticeSpacePartitioning(NULL, &ofs);

		ofs << "</body></html>"; ofs.close();
	}

	/*
	 * for i := 0 to n do begin
	 * 	 for j := 0 to n do begin
	 *     A[i][j] := A[i + 2][j] + A[i][j + 2];
	 *   end;
	 * end;
	 */
	void CRDMGraph_SpacePartition_PACTReview_AffineSet_Param() {
		test_RDMGraph_SetupGraph(1, 1, 1);

		cV[0]->setData(IP(4, IV_(4,  1,  0,  0,  0),
				             IV_(4, -1,  0,  1,  0),
				             IV_(4,  0,  1,  0,  0),
				             IV_(4,  0, -1,  1,  0)));
		dV[0]->setData(IPolyheder(2));

		g.addEdge(cV[0], dV[0], IAT(2, IV_(4,  1,  0,  0,  0),
				                       IV_(4,  0,  1,  0,  0)));
		g.addEdge(cV[0], dV[0], IAT(2, IV_(4,  1,  0,  0,  2),
				                       IV_(4,  0,  1,  0,  0)));
		g.addEdge(cV[0], dV[0], IAT(2, IV_(4,  1,  0,  0,  0),
				                       IV_(4,  0,  1,  0,  2)));

		ofstream ofs("examples/CRDMGraph_SpacePartition_PACTReview_AffineSet_Param.xhtml");
		ofs << "<html xmlns=\"http://www.w3.org/1999/xhtml\"><body>";

		performSpacePartitioning(NULL, &ofs);

/*		ofs << "<font size=\"+3\"><pre>"
		    << "for (i := 0 to 10) {\n"
		    << "  A[i + 1] = A[i - 2] * B[5*i];\n"
		    << "  B[5*i + 4] = A[4*i] + B[5*i + 4];\n"
		    << "  C[i + 3] = C[i + 7] + C[i + 3];\n"
		    << "}\n"
		    << "</pre></font>"*/;

		//encodeFWResultsInHtml(ofs, IntV(2, 0, 0), IntV(2, 10, 10));

		/*ASSERT(resultEdgesLattice[0][0]->getData() == ILR(1, 1, 1, 3, IV_(3, 1, 1, 0),
																 	  IV_(3, 0, 3, 0),
																	  IV_(3, 0, 0, 1)));
		ASSERT(resultEdgesLattice[0][1]->getData() == ILR(1, 1, 1, 3, IV_(3, 1, 0, 2),
						  											  IV_(3, 0, 1, 1),
						  											  IV_(3, 0, 0, 3)));
		ASSERT(resultEdgesLattice[0][1]->getData() == ILR(1, 1, 1, 3, IV_(3, 1, 0, 1),
						  											  IV_(3, 0, 1, 2),
						  											  IV_(3, 0, 0, 3)));
		ASSERT(resultEdgesLattice[1][1]->getData() == ILR(1, 1, 1, 3, IV_(3, 1, 1, 0),
																	  IV_(3, 0, 3, 0),
																	  IV_(3, 0, 0, 1)));*/

		//ofs << "<hr />";
		//determineDistanceLattices(&ofs);

		//distanceLattices[0]->print();

		//encodeDistanceLatticesInHtml(ofs, IntV(2, -7, -7), IntV(2,  7,  7));

		ofs << "</body></html>"; ofs.close();
	}

	/* Simple 1D uniform loopnest for lattice-based parallelism
	 * for ( i = 0 ; i <= n ; i++) {
	 *   A[i] := A[i] + B[i - 5];
	 *   B[i] := A[i - 1] * B[i];
	 * }
	 */
	void CRDMGraph_SpacePartition_1DUniformLattice() {
		test_RDMGraph_SetupGraph(2, 2, 1);

		cV[0]->setData(IPolyheder(2));
		dV[0]->setData(IPolyheder(1));
		cV[1]->setData(IPolyheder(2));
		dV[1]->setData(IPolyheder(1));

		g.addEdge(cV[0], dV[0], IAT(1, IV_(3, 1, 0,  0)));
		g.addEdge(cV[0], dV[1], IAT(1, IV_(3, 1, 0, -5)));
		g.addEdge(cV[1], dV[0], IAT(1, IV_(3, 1, 0, -1)));
		g.addEdge(cV[1], dV[1], IAT(1, IV_(3, 1, 0,  0)));

		ofstream ofs("examples/CRDMGraph_SpacePartition_1DUniformLattice.xhtml");
		ofs << "<html xmlns=\"http://www.w3.org/1999/xhtml\"><body>";

		performLatticeSpacePartitioning(NULL, &ofs);

		//encodeFWResultsInHtml(ofs, IntV(2, 0, 0), IntV(2, 15, 15));

		ASSERT(resultEdgesLattice[0][0]->getData() == ILR(1, 1, 2, 4, IV_(4, 1, 1, 0, 0),
																	 		  IV_(4, 0, 6, 0, 0),
																			  IV_(4, 0, 0, 1, 0),
																			  IV_(4, 0, 0, 0, 1)));
		ASSERT(resultEdgesLattice[0][1]->getData() == ILR(1, 1, 2, 4, IV_(4, 1, 0, 0, 5),
						  IV_(4, 0, 1, 0, 1),
						  IV_(4, 0, 0, 1, 0),
						  										 IV_(4, 0, 0, 0, 6)));
		ASSERT(resultEdgesLattice[1][0]->getData() == ILR(1, 1, 2, 4, IV_(4, 1, 0, 0, 1),
																	 IV_(4, 0, 1, 0, 5),
																	 IV_(4, 0, 0, 1, 0),
																	 IV_(4, 0, 0, 0, 6)));
		ASSERT(resultEdgesLattice[1][1]->getData() == ILR(1, 1, 2, 4, IV_(4, 1, 1, 0, 0),
																	 IV_(4, 0, 6, 0, 0),
																	 IV_(4, 0, 0, 1, 0),
																	 IV_(4, 0, 0, 0, 1)));

		ofs << "</body></html>"; ofs.close();
	}

	/* Simple 1D uniform loopnest for lattice-based parallelism
	 * for ( i = 0 ; i <= n ; i++) {
	 *   A[i] := A[i] + B[i - 5] + C[i + 3];
	 *   B[i] := A[i - 1] * B[i];
	 * 	 C[i] := A[i - 3] + B[i + 1];
	 * }
	 */
	void CRDMGraph_SpacePartition_1DUniformLattice_3Statements() {
		test_RDMGraph_SetupGraph(3, 3, 1);

		cV[0]->setData(IPolyheder(2));
		cV[1]->setData(IPolyheder(2));
		cV[2]->setData(IPolyheder(2));
		dV[0]->setData(IPolyheder(2));
		dV[1]->setData(IPolyheder(2));
		dV[2]->setData(IPolyheder(2));

		g.addEdge(cV[0], dV[0], IAT(1, IV_(3, 1, 0,  0)));
		//g.addEdge(cV[0], dV[0], iat_im6);
		g.addEdge(cV[0], dV[1], IAT(1, IV_(3, 1, 0, -5)));
		g.addEdge(cV[1], dV[0], IAT(1, IV_(3, 1, 0, -1)));
		g.addEdge(cV[1], dV[1], IAT(1, IV_(3, 1, 0,  0)));
		g.addEdge(cV[2], dV[0], IAT(1, IV_(3, 1, 0, -3)));
		g.addEdge(cV[2], dV[1], IAT(1, IV_(3, 1, 0,  1)));
		g.addEdge(cV[0], dV[2], IAT(1, IV_(3, 1, 0,  3)));

		performLatticeSpacePartitioning();

		ofstream ofs("CRDMGraph_SpacePartition_1DUniformLattice_3Statements.xhtml");
		ofs << "<html xmlns=\"http://www.w3.org/1999/xhtml\"><body>";

		encodeFWResultsInHtml(ofs, IntV(2, 0, 0), IntV(2, 11, 11));

		ofs << "<hr />";

		determineDistanceLattices();

		encodeDistanceLatticesInHtml(ofs, IntV(1, -7), IntV(1, 7));

		ofs << "</body></html>";
		/*ASSERT((resultEdges[0][1].size() == 1) && (resultEdges[0][1][0]->getData() == imr_imj_imj_m1));
		ASSERT((resultEdges[1][0].size() == 1) && (resultEdges[1][0][0]->getData() == imr_imj_imj_1));
		ASSERT((resultEdges[0][0].size() == 1) && (resultEdges[0][0][0]->getData() == imr_imj_imj));
		ASSERT((resultEdges[1][1].size() == 1) && (resultEdges[1][1][0]->getData() == imr_imj_imj));*/
	}

	void CRDMGraph_AffineSet_SpacePartition_2DUniformLattice_YuRot() {
		test_RDMGraph_SetupGraph(2, 2, 0);

		cV[0]->setData(IP(4, IV_(3,  1,  0,  9),
	                         IV_(3, -1,  0, -9),
	                         IV_(3,  0,  1,  8),
	                         IV_(3,  0, -1, -8)));
		cV[1]->setData(IP(4, IV_(3,  1,  0,  9),
	                         IV_(3, -1,  0, -9),
	                         IV_(3,  0,  1,  8),
	                         IV_(3,  0, -1, -8)));
		dV[0]->setData(IPolyheder(2));
		dV[1]->setData(IPolyheder(2));

		g.addEdge(cV[0], dV[0], IAT(2, IV_(3,  4, -1,  3),
				                       IV_(3,  2,  1, -2)));
		g.addEdge(cV[0], dV[0], IAT(2, IV_(3,  1,  1, -1),
				                       IV_(3,  1, -1,  2)));
		g.addEdge(cV[1], dV[0], IAT(2, IV_(3,  1,  0, 2),
				                       IV_(3,  1,  0, -2)));
		g.addEdge(cV[1], dV[1], IAT(2, IV_(3,  1,  0,  0),
				                       IV_(3,  0,  1,  0)));
		g.addEdge(cV[1], dV[1], IAT(2, IV_(3,  1,  0,  0),
                					   IV_(3,  0,  1, -2)));

		ofstream ofs("examples/CRDMGraph_AffineSet_SpacePartition_2DUniformLattice_YuRot.xhtml");
		ofs << "<html xmlns=\"http://www.w3.org/1999/xhtml\"><body>";

		performSpacePartitioning(NULL, &ofs);

		ofs << "</body></html>"; ofs.close();

		/*ASSERT((resultEdges[0][1].size() == 1) && (resultEdges[0][1][0]->getData() == imr_imj_imj_m1));
		ASSERT((resultEdges[1][0].size() == 1) && (resultEdges[1][0][0]->getData() == imr_imj_imj_1));
		ASSERT((resultEdges[0][0].size() == 1) && (resultEdges[0][0][0]->getData() == imr_imj_imj));
		ASSERT((resultEdges[1][1].size() == 1) && (resultEdges[1][1][0]->getData() == imr_imj_imj));*/
	}

	void CRDMGraph_AffineLattice_SpacePartition_2DUniformLattice_YuRot() {
		test_RDMGraph_SetupGraph(2, 2, 0);

		cV[0]->setData(IP(4, IV_(3,  1,  0,  9),
	                         IV_(3, -1,  0, -9),
	                         IV_(3,  0,  1,  8),
	                         IV_(3,  0, -1, -8)));
		cV[1]->setData(IP(4, IV_(3,  1,  0,  9),
	                         IV_(3, -1,  0, -9),
	                         IV_(3,  0,  1,  8),
	                         IV_(3,  0, -1, -8)));
		dV[0]->setData(IPolyheder(2));
		dV[1]->setData(IPolyheder(2));

		g.addEdge(cV[0], dV[0], IAT(2, IV_(3,  4, -1,  3),
				                       IV_(3,  2,  1, -2)));
		g.addEdge(cV[0], dV[0], IAT(2, IV_(3,  1,  1, -1),
				                       IV_(3,  1, -1,  2)));
		g.addEdge(cV[1], dV[0], IAT(2, IV_(3,  1,  0, 2),
				                       IV_(3,  1,  0, -2)));
		g.addEdge(cV[1], dV[1], IAT(2, IV_(3,  1,  0,  0),
				                       IV_(3,  0,  1,  0)));
		g.addEdge(cV[1], dV[1], IAT(2, IV_(3,  1,  0,  0),
                					   IV_(3,  0,  1, -2)));

		ofstream ofs("examples/CRDMGraph_AffineLattice_SpacePartition_2DUniformLattice_YuRot.xhtml");
		ofs << "<html xmlns=\"http://www.w3.org/1999/xhtml\"><body>";

		performLatticeSpacePartitioning(NULL, &ofs);

		ofs << "</body></html>"; ofs.close();

		/*ASSERT((resultEdges[0][1].size() == 1) && (resultEdges[0][1][0]->getData() == imr_imj_imj_m1));
		ASSERT((resultEdges[1][0].size() == 1) && (resultEdges[1][0][0]->getData() == imr_imj_imj_1));
		ASSERT((resultEdges[0][0].size() == 1) && (resultEdges[0][0][0]->getData() == imr_imj_imj));
		ASSERT((resultEdges[1][1].size() == 1) && (resultEdges[1][1][0]->getData() == imr_imj_imj));*/
	}

	/* From: Partitioning Loops with Variable Dependence Distance (Yu, D'Hollander)
	 */
	void CRDMGraph_SpacePartition_Yu_1() {
		test_RDMGraph_SetupGraph(1, 1, 0);

		cV[0]->setData(IPolyheder(2));
		dV[0]->setData(IPolyheder(2));

		g.addEdge(cV[0], dV[0], IAT(2, IV_(3, 3, 0,  1),
				                       IV_(3, 2, 1, -1)));
		g.addEdge(cV[0], dV[0], IAT(2, IV_(3, 1, 0,  3),
				                       IV_(3, 0, 1,  1)));

		performLatticeSpacePartitioning();

		determineDistanceLattices();

		ofstream ofs("CRDMGraph_SpacePartition_Yu_1.xhtml");
		ofs << "<html xmlns=\"http://www.w3.org/1999/xhtml\"><body>";

		encodeDistanceLatticesInHtml(ofs, IntV(2, -7, -7), IntV(2, 7, 7));

		/*ASSERT((resultEdgesLattice[0][0].size() == 1)
				&& (resultEdgesLattice[0][0][0]->getData() == ILR(1, 1, 1, 3, IV_(3, 1, 1, 0),
																			  IV_(3, 0, 6, 0),
																			  IV_(3, 0, 0, 1))));

		ofs << "</table>";*/

		ofs << "</body></html>"; ofs.close();
	}

	/* From: Partitioning Loops with Variable Dependence Distance (Yu, D'Hollander)
	 */
	void CRDMGraph_SpacePartition_Yu_2() {
		test_RDMGraph_SetupGraph(1, 1, 0);

		cV[0]->setData(IPolyheder(2));
		dV[0]->setData(IPolyheder(2));

		g.addEdge(cV[0], dV[0], IAT(2, IV_(3,  4, -1,  3),
				                       IV_(3,  2,  1, -2)));
		g.addEdge(cV[0], dV[0], IAT(2, IV_(3,  1,  1, -1),
				                       IV_(3,  1, -1,  2)));

		performLatticeSpacePartitioning();

		determineDistanceLattices();

		ofstream ofs("examples/CRDMGraph_SpacePartition_Yu_2.xhtml");
		ofs << "<html xmlns=\"http://www.w3.org/1999/xhtml\"><body>";

		encodeDistanceLatticesInHtml(ofs, IntV(2, -7, -7), IntV(2, 7, 7));

		/*ASSERT((resultEdgesLattice[0][0].size() == 1)
				&& (resultEdgesLattice[0][0][0]->getData() == ILR(1, 1, 1, 3, IV_(3, 1, 1, 0),
																			  IV_(3, 0, 6, 0),
																			  IV_(3, 0, 0, 1))));

		ofs << "</table>";*/

		ofs << "</body></html>"; ofs.close();
	}

	/**
	 * Test Lim & Lam (Communication-Free Parallelisation via Affine Transformations)
	 * - To show that some sets of dependences lead to processor spaces which are not full-dimensional
	 * for ( i = 0 ; i <= n ; i++) {
	 * 	 for ( j = 0 ; j <= n ; j++) {
	 * 		A[i][j][1] = A[i][j][1] + ...
	 *   }
	 * }
	 * for ( i = 0 ; i <= n ; i++) {
	 * 	 for ( j = 0 ; j <= n ; j++) {
	 * 		A[i][1][j] = A[i][1][j] + ...
	 *   }
	 * }
	 * for ( i = 0 ; i <= n ; i++) {
	 * 	 for ( j = 0 ; j <= n ; j++) {
	 * 		A[1][i][j] = A[1][i][j] + ...
	 *   }
	 * }
	 *
	 * Dependencies:
	 * - SpaceVector intersection
	 */
	void CRDMGraph_AffineSet_SpacePartition_LimLam_NonFullDimensional() {
		test_RDMGraph_SetupGraph(3, 1, 1);

		cV[0]->setData(IP(4, IV_(4,  1,  0,  0,  0),
	             		     IV_(4, -1,  0,  1,  0),
	             		     IV_(4,  0,  1,  0,  0),
	             		     IV_(4,  0, -1,  1,  0)));

		cV[1]->setData(IP(4, IV_(4,  1,  0,  0,  0),
	             		     IV_(4, -1,  0,  1,  0),
	             		     IV_(4,  0,  1,  0,  0),
	             		     IV_(4,  0, -1,  1,  0)));

		cV[2]->setData(IP(4, IV_(4,  1,  0,  0,  0),
	             		     IV_(4, -1,  0,  1,  0),
	             		     IV_(4,  0,  1,  0,  0),
	             		     IV_(4,  0, -1,  1,  0)));
		dV[0]->setData(IPolyheder(3));

		g.addEdge(cV[0], dV[0], IAT(3, IV_(4,  1,  0,  0,  0),
				                       IV_(4,  0,  1,  0,  0),
				                       IV_(4,  0,  0,  0,  1)));

		g.addEdge(cV[1], dV[0], IAT(3, IV_(4,  1,  0,  0,  0),
				                       IV_(4,  0,  0,  0,  1),
				                       IV_(4,  0,  1,  0,  0)));

		g.addEdge(cV[2], dV[0], IAT(3, IV_(4,  0,  0,  0,  1),
				                       IV_(4,  1,  0,  0,  0),
				                       IV_(4,  0,  1,  0,  0)));

		ofstream ofs("examples/CRDMGraph_AffineSet_SpacePartition_LimLam_NonFullDimensional.xhtml");
		ofs << "<html xmlns=\"http://www.w3.org/1999/xhtml\"><body>";

		performSpacePartitioning(NULL, &ofs);

/*
		ASSERT((resultEdges[0][0].size() == 1) && (resultEdges[0][0][0]->getData() == imr_i_i_0__j_j_0));
		ASSERT((resultEdges[1][1].size() == 1) && (resultEdges[1][1][0]->getData() == imr_i_i_0__j_j_0));
		ASSERT((resultEdges[2][2].size() == 1) && (resultEdges[2][2][0]->getData() == imr_i_i_0__j_j_0));

		ASSERT((resultEdges[0][1].size() == 1) && (resultEdges[0][1][0]->getData() == imr_i_i_0__j_0_m1___0_j_1));
		ASSERT((resultEdges[1][0].size() == 1) && (resultEdges[1][0][0]->getData() == imr_i_i_0__j_0_m1___0_j_1));

		ASSERT((resultEdges[0][2].size() == 1) && (resultEdges[0][2][0]->getData() == imr_j_i_0__i_0_m1___0_j_1));
		ASSERT((resultEdges[2][0].size() == 1) && (resultEdges[2][0][0]->getData() == imr_i_j_0__j_0_m1___0_i_1));

		ASSERT((resultEdges[1][2].size() == 1) && (resultEdges[1][2][0]->getData() == imr_j_j_0__i_0_m1___0_i_1));
		ASSERT((resultEdges[2][1].size() == 1) && (resultEdges[2][1][0]->getData() == imr_j_j_0__i_0_m1___0_i_1));
*/
		ofs << "</body></html>"; ofs.close();
	}

	/**
	 * Test Lim & Lam (Communication-Free Parallelisation via Affine Transformations)
	 * - To show that some sets of dependences lead to processor spaces which are not full-dimensional
	 * for ( i = 0 ; i <= n ; i++) {
	 * 	 for ( j = 0 ; j <= n ; j++) {
	 * 		A[i][j][1] = A[i][j][1] + ...
	 *   }
	 * }
	 * for ( i = 0 ; i <= n ; i++) {
	 * 	 for ( j = 0 ; j <= n ; j++) {
	 * 		A[i][1][j] = A[i][1][j] + ...
	 *   }
	 * }
	 * for ( i = 0 ; i <= n ; i++) {
	 * 	 for ( j = 0 ; j <= n ; j++) {
	 * 		A[1][i][j] = A[1][i][j] + ...
	 *   }
	 * }
	 *
	 * Dependencies:
	 * - SpaceVector intersection
	 */
	void CRDMGraph_AffineLattice_SpacePartition_LimLam_NonFullDimensional() {
		test_RDMGraph_SetupGraph(3, 1, 1);

		cV[0]->setData(IP(4, IV_(4,  1,  0,  0,  0),
	             		     IV_(4, -1,  0,  1,  0),
	             		     IV_(4,  0,  1,  0,  0),
	             		     IV_(4,  0, -1,  1,  0)));

		cV[1]->setData(IP(4, IV_(4,  1,  0,  0,  0),
	             		     IV_(4, -1,  0,  1,  0),
	             		     IV_(4,  0,  1,  0,  0),
	             		     IV_(4,  0, -1,  1,  0)));

		cV[2]->setData(IP(4, IV_(4,  1,  0,  0,  0),
	             		     IV_(4, -1,  0,  1,  0),
	             		     IV_(4,  0,  1,  0,  0),
	             		     IV_(4,  0, -1,  1,  0)));
		dV[0]->setData(IPolyheder(3));

		g.addEdge(cV[0], dV[0], IAT(3, IV_(4,  1,  0,  0,  0),
				                       IV_(4,  0,  1,  0,  0),
				                       IV_(4,  0,  0,  0,  1)));

		g.addEdge(cV[1], dV[0], IAT(3, IV_(4,  1,  0,  0,  0),
				                       IV_(4,  0,  0,  0,  1),
				                       IV_(4,  0,  1,  0,  0)));

		g.addEdge(cV[2], dV[0], IAT(3, IV_(4,  0,  0,  0,  1),
				                       IV_(4,  1,  0,  0,  0),
				                       IV_(4,  0,  1,  0,  0)));

		ofstream ofs("examples/CRDMGraph_AffineLattice_SpacePartition_LimLam_NonFullDimensional.xhtml");
		ofs << "<html xmlns=\"http://www.w3.org/1999/xhtml\"><body>";

		performLatticeSpacePartitioning(NULL, &ofs);

		ofs << "</body></html>"; ofs.close();
	}
}

cute::suite* Test_RDMGraph_runSuite(){
	cute::suite &s = *(new cute::suite("RDMGraph"));

	/*s.push_back(CUTE(Test_RDMGraph::CRDMGraph_SpacePartition_Trivial2DParallelism));

	s.push_back(CUTE(Test_RDMGraph::CRDMGraph_AffineSet_SpacePartition_LimLam));
	s.push_back(CUTE(Test_RDMGraph::CRDMGraph_AffineLattice_SpacePartition_LimLam));

	s.push_back(CUTE(Test_RDMGraph::CRDMGraph_AffineSet_SpacePartition_LimLam_NonFullDimensional));
	s.push_back(CUTE(Test_RDMGraph::CRDMGraph_AffineLattice_SpacePartition_LimLam_NonFullDimensional));

	s.push_back(CUTE(Test_RDMGraph::CRDMGraph_SpacePartition_PACTReview_AffineLattice_Param));
	s.push_back(CUTE(Test_RDMGraph::CRDMGraph_SpacePartition_PACTReview_AffineSet_Param));
	s.push_back(CUTE(Test_RDMGraph::CRDMGraph_SpacePartition_1DUniformLattice));
	s.push_back(CUTE(Test_RDMGraph::CRDMGraph_AffineSet_SpacePartition_1DUniformLattice_NonUniform));
	s.push_back(CUTE(Test_RDMGraph::CRDMGraph_AffineLattice_SpacePartition_1DUniformLattice_NonUniform));
	s.push_back(CUTE(Test_RDMGraph::CRDMGraph_SpacePartition_1DUniformLattice_NonUniform3));
	s.push_back(CUTE(Test_RDMGraph::CRDMGraph_SpacePartition_1DUniformLattice_NonUniform3_AffineSet));
	s.push_back(CUTE(Test_RDMGraph::CRDMGraph_AffineSet_SpacePartition_2DUniformLattice_YuRot));
	s.push_back(CUTE(Test_RDMGraph::CRDMGraph_AffineLattice_SpacePartition_2DUniformLattice_YuRot));*/
	/*s.push_back(CUTE(Test_RDMGraph::CRDMGraph_SpacePartition_1DUniformLattice_SingleStatement));
	s.push_back(CUTE(Test_RDMGraph::CRDMGraph_SpacePartition_1DUniformLattice_SingleStatement_NonUniform));
	s.push_back(CUTE(Test_RDMGraph::CRDMGraph_SpacePartition_1DUniformLattice_NonUniform));
	s.push_back(CUTE(Test_RDMGraph::CRDMGraph_SpacePartition_1DUniformLattice_NonUniform2));
	s.push_back(CUTE(Test_RDMGraph::CRDMGraph_SpacePartition_1DUniformLattice_NonUniform3));
	s.push_back(CUTE(Test_RDMGraph::CRDMGraph_SpacePartition_1DUniformLattice_NonUniform_3Statement));
	s.push_back(CUTE(Test_RDMGraph::CRDMGraph_SpacePartition_1DUniformLattice_3Statements));
	s.push_back(CUTE(Test_RDMGraph::CRDMGraph_SpacePartition_2DUniformLattice_YuRot));
	s.push_back(CUTE(Test_RDMGraph::CRDMGraph_SpacePartition_Yu_1));
	s.push_back(CUTE(Test_RDMGraph::CRDMGraph_SpacePartition_Yu_2));*/

	//s.push_back(CUTE(Test_RDMGraph::CRDMGraph_SpacePartition_2DUniformLattice));

	return &s;
}

#endif
