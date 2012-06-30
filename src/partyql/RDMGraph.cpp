/* OBSOLETE */

#ifndef _RDMGRAPH_CPP
#define _RDMGRAPH_CPP

#include "RDMGraph.h"

#include "../utils/ToString.h"
#include "../utils/SVG.h"
#include "../utils/SVGELattice.h"
#include "basicmath/scalar/Integer.h"
#include "basicmath/SetRelation.h"

namespace AlgoTrans {
	template <class R>
	typename CModuleRADGraph<R>::T CRDMGraph<R>::calculateAffinelyApproximatedRSDG() {
		typedef typename CModuleRelation<R>::T ModuleRelation;
		typedef typename CModuleRADGraph<R>::T RADGraph;
		typedef typename CModuleRADGVertex<R>::T RADGVertex;

		/*for (unsigned int v = 0; v < computationVertices.size(); v++) { // We need to implement the hierarchical parameter-mapping before we can improve this
			ASSERT((v == 0) || (parmCount == computationVertices[v]->getData().getParameterCount()));
			parmCount = computationVertices[v]->getData().getParameterCount();
		}*/

		if (this->osHtml != NULL) {
			*this->osHtml << "Parameters:" << parmCount << "<br/><b>Basic edge relations (each represented as vector space basis):</b><br/><table>"
			<< "<tr><td>Variable</td><td>Left statement</td><td>Right statement</td><td>Pairwise relation</td><td>Symmetric relation</td></tr>";
		}
		// Calculate vertex-pair spaces (result stored in new graph)
		RADGraph radg = RADGraph();
		radg.osHtml = this->osHtml;
		for (unsigned int v = 0; v < computationVertices.size(); v++) radg.addVertex(computationVertices[v]->getData());
		for (unsigned int v = 0; v < computationVertices.size(); v++) { CRDMGVertex& vxV = *computationVertices[v];
			RADGVertex& resultVxV = radg(v);
			for (unsigned int w = v; w < computationVertices.size(); w++) { CRDMGVertex& vxW = *computationVertices[w];
				RADGVertex& resultVxW = radg(w);

				int vxVUDim = vxV.getData().getSpaceDimension() - parmCount;
				int vxWUDim = vxW.getData().getSpaceDimension() - parmCount;

				//int localUnknownSpaceDim = vxV.getData().getSpaceDimension() + vxW.getData().getSpaceDimension() - 2*parmCount;
				//int localSpaceDim = localUnknownSpaceDim + parmCount + 1;

				// Calculate solution-space for a partition for vertex-pair (v, w)
				vector<CModuleRelation<R> > accessPairModules = vector<CModuleRelation<R> >();
				for (unsigned int a = 0; a < dataVertices.size(); a++) { CRDMGVertex& vxData = *dataVertices[a];
					const vector<CRDMGEdge*> edgesVToData = getEdgesBetween(vxV, vxData);
					const vector<CRDMGEdge*> edgesWToData = getEdgesBetween(vxW, vxData);
					for (unsigned int e = 0; e < edgesVToData.size(); e++) {
						CAffineTransformation<R>& atE = edgesVToData[e]->getData();
						ASSERT(atE.getDestinationSpaceDimension() == vxData.getData().getSpaceDimension());

						for (unsigned int f = 0; f < edgesWToData.size(); f++) if ((v != w) || (e <= f)) {
							CAffineTransformation<R>& atF = edgesWToData[f]->getData();
							//ASSERT(&atE != &atF);

							bool osHtmlAct = (this->osHtml != NULL);
							if (osHtmlAct) *this->osHtml << "<tr><td>" << a << "</td><td>" << v << "</td><td>" << w << "</td>";

							// Fill in null-space vectors obtained through data-access-equation
							ModuleRelation mLR = ModuleRelation::fromConstraints(vxVUDim, vxWUDim, parmCount + 1,
								atE.getSubMatrixColumns(0, vxVUDim) << -atF.getSubMatrixColumns(0, vxWUDim)
								<< (atE.getSubMatrixColumns(vxVUDim, parmCount + 1) - atF.getSubMatrixColumns(vxWUDim, parmCount + 1))
							);
							ModuleRelation mRL = ModuleRelation::fromConstraints(vxVUDim, vxWUDim, parmCount + 1,
								atF.getSubMatrixColumns(0, vxVUDim) << -atE.getSubMatrixColumns(0, vxWUDim)
								<< (atF.getSubMatrixColumns(vxVUDim, parmCount + 1) - atE.getSubMatrixColumns(vxWUDim, parmCount + 1))
							);
							radg.addEdge(resultVxV, resultVxW, mLR);
							radg.addEdge(resultVxW, resultVxV, mRL);

							if (osHtmlAct) *this->osHtml << "<td>" << ModuleRelation(mLR).toStringHtml() << "</td>";
							if (osHtmlAct) *this->osHtml << "<td>" << ModuleRelation(mRL).toStringHtml() << "</td></tr>";
						}
					}
				}

				/*if (v == w) { // Intersect accessPairModules with coefficient equality constraints
					// THINK: So how can this be done more efficiently?
					// THINK: And what about the constants?
					// THINK: Should the unknown dimensions be treated the same as the constants above for v = w?
					CModuleRelation<R> m = CModuleRelation<R>(vxVUDim, vxVUDim, parmCount + 1);
					for (int d = 0; d < vxV.getData().getSpaceDimension() - parmCount; d++) {
						m.addGeneratingVertex(new CVector<R>(
							CVector<R>::getUnitVector(vxVUDim, d) << CVector<R>::getUnitVector(vxVUDim, d) << CVector<R>::getZeroVector(parmCount + 1)
						));
					}
					radg.addEdge(resultVxV, resultVxW, m);
					radg.addEdge(resultVxW, resultVxV, m);
				}*/
			}
		}

		if (this->osHtml != NULL) *this->osHtml << "</table>";

		return radg;
	}

	template <class R>
	typename CLatticeRADGraph<R>::T CRDMGraph<R>::calculateFlatLatticeApproximatedRSDG() {
		typedef typename CLatticeRelation<R>::T LatticeRelation;
		typedef typename CLatticeRADGraph<R>::T RADGraph;
		typedef typename CLatticeRADGVertex<R>::T RADGVertex;

		/*for (unsigned int v = 0; v < computationVertices.size(); v++) { // We need to implement the hierarchical parameter-mapping before we can improve this
			ASSERT((v == 0) || (parmCount == computationVertices[v]->getData().getParameterCount()));
			parmCount = computationVertices[v]->getData().getParameterCount();
		}*/

		vector<string>* dimNamesI = NULL;
		vector<string>* dimNamesIp = NULL;
		vector<string>* dimNamesIIp = NULL;
		if (this->osLatex != NULL) {
			*this->osLatex << "\\begin{tabular}{| l | c | c |}\n";
			*this->osLatex << "\\hline\n";
			*this->osLatex
			<< "$\\aVariableAccessFunctionA(p)$, $\\aVariableAccessFunctionB(q)$ "
			<< "& $\\aVariableAccessFunctionA(p) - \\aVariableAccessFunctionB(q) \\eq 0$ "
			<< "& Solutions " << "\\\\\n";
		}
		if (this->osHtml != NULL) {
			*this->osHtml << "Parameters:" << parmCount << "<br/><b>Basic edge relations:</b><br/><table>"
			<< "<tr><td>Variable</td><td>Left statement</td><td>Right statement</td><td>Relation equation</td><td>Symmetric relation equation</td></tr>";
		}

		/*if (this->osLatex != NULL) {
			dimNamesI = new vector<string>();
			dimNamesI->push_back("p");

			dimNamesIp = new vector<string>();
			dimNamesIp->push_back("q");

			dimNamesIIp = new vector<string>();
			dimNamesIIp->push_back("p");
			dimNamesIIp->push_back("q");
		}*/
		int iTitle = 0;
		string abcdef = "abcdef";
		unsigned int lV = 0, lW = 0, lA = 0;
		bool valid = false;
		int svgIx = 0;

		// Calculate vertex-pair spaces (result stored in new graph)
		RADGraph radg = RADGraph();
		radg.osHtml = this->osHtml;
		for (unsigned int v = 0; v < computationVertices.size(); v++) radg.addVertex(computationVertices[v]->getData());
		for (unsigned int v = 0; v < computationVertices.size(); v++) { CRDMGVertex& vxV = *computationVertices[v];
			RADGVertex& resultVxV = radg(v);
			for (unsigned int w = v; w < computationVertices.size(); w++) { CRDMGVertex& vxW = *computationVertices[w];
				RADGVertex& resultVxW = radg(w);

				int vxVUDim = vxV.getData().getSpaceDimension() - parmCount;
				int vxWUDim = vxW.getData().getSpaceDimension() - parmCount;

				// Calculate dependence lattice induced by the common memory accesses of vertex-pair (v, w)
				vector<CModuleRelation<R> > accessPairModules = vector<CModuleRelation<R> >();
				for (unsigned int a = 0; a < dataVertices.size(); a++) { CRDMGVertex& vxData = *dataVertices[a];
					const vector<CRDMGEdge*> edgesVToData = getEdgesBetween(vxV, vxData);
					const vector<CRDMGEdge*> edgesWToData = getEdgesBetween(vxW, vxData);
					bool firstF = true;
					for (unsigned int e = 0; e < edgesVToData.size(); e++) {
						CAffineTransformation<R>& atE = edgesVToData[e]->getData();
						ASSERT(atE.getDestinationSpaceDimension() == vxData.getData().getSpaceDimension());

						for (unsigned int f = 0; f < edgesWToData.size(); f++) if ((v != w) || (e <= f)) {
							CAffineTransformation<R>& atF = edgesWToData[f]->getData();

							bool osLatexAct = (this->osLatex != NULL) && (&atE != &atF);
							bool osHtmlAct = (this->osHtml != NULL);
							if (osHtmlAct) *this->osHtml << "<tr><td>" << a << "</td><td>" << v << "</td><td>" << w << "</td>";
							if (osLatexAct) {
								*this->osLatex << "\\hline\n";
								if (!valid || (lV != v) || (lW != w) || (lA != a)) {
									*this->osLatex << "\\multicolumn{3}{|l|}{" <<
									"Array " << (a == 0 ? "A" : "B") << " accessed by statement pair (" << v << ", " << w << "):";
									//"Statement pair (" << v << ", " << w << ") accessing array " << a << ":}\\\\";
									*this->osLatex << "} \\\\ \\hline\n";
								}
								lV = v; lW = w; lA = a;
								valid = true;

								*this->osLatex << "$\\begin{array}{ll}\nf(p) \\eq " << atE.toStringLatex(dimNamesI) << "\\\\";
								*this->osLatex << "g(q) \\eq " << atF.toStringLatex(dimNamesIp) << "\n\\end{array}$\n";
							}
							firstF = false;

							// Fill in null-space vectors obtained through data-access-equation
							CMatrix<R> mLR = CMatrix<R>(vxVUDim + vxWUDim + parmCount + 1);
							CMatrix<R> mRL = CMatrix<R>(vxVUDim + vxWUDim + parmCount + 1);
							for (int d = 0; d < vxData.getData().getSpaceDimension(); d++) {
								// XXX: We may want to avoid unused dimensions for computational efficiency
								mLR.addRow(
									atE[d].getSubVector(0, vxVUDim) << -atF[d].getSubVector(0, vxWUDim)
									<< (atE[d].getSubVector(vxVUDim, parmCount + 1) - atF[d].getSubVector(vxWUDim, parmCount + 1)));
								mRL.addRow(
									atF[d].getSubVector(0, vxWUDim) << -atE[d].getSubVector(0, vxVUDim)
									<< (atF[d].getSubVector(vxWUDim, parmCount + 1) - atE[d].getSubVector(vxVUDim, parmCount + 1)));
							}

							if (osLatexAct) *this->osLatex << " & $" << CAffineTransformation<R>(mLR).toStringLatex(dimNamesIIp) << " \\eq 0$";

							if (osHtmlAct) *this->osHtml << "<td>" << CAffineTransformation<R>(mLR).toStringHtml(dimNamesIIp) << "</td>";
							if (osHtmlAct) *this->osHtml << "<td>" << CAffineTransformation<R>(mRL).toStringHtml(dimNamesIIp) << "</td></tr>";

							// HNFing here for debugging convenience -- not necessary (performance?)
							// mLR can most likely be easily obtained from mRL without double HNFing or !ing
							LatticeRelation gLR = LatticeRelation::fromConstraints(vxVUDim, vxWUDim, parmCount + 1, mLR);
							LatticeRelation gRL = LatticeRelation::fromConstraints(vxWUDim, vxVUDim, parmCount + 1, mRL);

							if (osLatexAct) {
								CVector<R> exemplar = gLR.getAffineExemplar();
								CMatrix<R> linLat = gLR.getLinearLattice().getMatrix().getTransposedMatrix();

 								CSVGELattice<R>* lat = new CSVGELattice<R>(gLR, latticeLowerBounds, latticeUpperBounds, latticeBasicWidth, true);
 								string pq = "pq";
 								lat->horAxisTitle = "p";
 								lat->verAxisTitle = "q";
 								string title = "-"; title[0] = abcdef[iTitle++];
 								lat->title = "(" + title + ")";
								int svgID = -1;
								for (int q = 0; q < svgIx; q++) {
									if (gLR == *(resultingLattices[q])) {
										svgID = q;
										iTitle--;
										delete lat;
										break;
									}
								}
								if (svgID == -1) {
									resultingSVGLattices.push_back(lat);
									resultingLattices.push_back(new CLattice<R>(CLattice<R>::fromGenerators(gLR.getMatrix())));
									svgID = svgIx++;
								}

								string abcde = "abcde";
								*this->osLatex << " & Figure \\ref{figureInitialLattices} (" << abcde[svgID] << ") \\\\";
								//*this->osLatex << " & $\\Lattice{" << linLat.toStringLatex() << ", " << exemplar.toStringLatex() << "}$ \\\\";
							}

							radg.addEdge(resultVxV, resultVxW, gLR);
							radg.addEdge(resultVxW, resultVxV, gRL);
						}
					}
				}

				/*if ((v == w)) { // Does this improve performance?
					CLatticeRelation<R> m = CLatticeRelation<R>(vxVUDim, vxVUDim, parmCount + 1);
					for (int d = 0; d < vxV.getData().getSpaceDimension() - parmCount; d++) {
						m.addGeneratingVertex(new CVector<R>(
							CVector<R>::getUnitVector(vxVUDim, d)
							<< CVector<R>::getUnitVector(vxVUDim, d)
							<< CVector<R>::getZeroVector(parmCount + 1)
						));
					}
					radg.addEdge(resultVxV, resultVxW, m);
					radg.addEdge(resultVxW, resultVxV, m);
				}*/
			}
		}

		if (this->osHtml != NULL) *this->osHtml << "</table>";

		if (this->osLatex != NULL) *this->osLatex << "\\hline\n";
		if (this->osLatex != NULL) *this->osLatex << "\\end{tabular}\n";

		return radg;
	}

	template <class R>
	CompGraphSpacePartitioner<R>* CRDMGraph<R>::performLatticeBasedSpacePartitioning() {
		CompGraphSpacePartitioner<R>* sp = new CompGraphSpacePartitioner<R>(this);

		sp.performLatticeSpacePartitioning();

		return sp;
	}

	template <class R>
	string CRDMGraph<R>::toStringHtml() {
		string result = "<table><tr><td>Statement</td><td>Variable</td><td>Reference</td></tr>";
		for (int p = 0; p < getComputationVertexCount(); p++) { CRDMGVertex& vxP = getComputationVertex(p);
			for (int d = 0; d < getDataVertexCount(); d++) { CRDMGVertex& vxD = getDataVertex(d);
				vector<CRDMGEdge* > edges = getEdgesBetween(vxP, vxD);
				for (unsigned int e = 0; e < edges.size(); e++) { CRDMGEdge& egE = *edges[e];
					result += "<tr><td>" + CInteger(p).toString() + "</td><td>" + CInteger(d).toString() + "</td><td>";
					result += egE.getData().toStringHtml();
					result += "</td></tr>";
				}
			}
		}
		result += "</table>";

		result += "<table><tr><td>Statement</td><td>Iteration Space Constraint Matrix</td></tr>";
		for (int p = 0; p < getComputationVertexCount(); p++) { CRDMGVertex& vxP = getComputationVertex(p);
			result += "<tr><td>" + CInteger(p).toString() + "</td><td>";
			result += vxP.getData().toStringHtml();
			result += "</td></tr>";
		}
		result += "</table>";

		return result;
	}
}

#endif /* _RDMGRAPH_CPP */
