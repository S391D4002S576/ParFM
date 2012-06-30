#ifndef COMPUTATION_H_
#define COMPUTATION_H_

#include <vector>
using std::vector;

#include <sstream>
using std::ostream;

#include <string>
using std::string;


#include "Declarations.h"

//#include "basicmath/PolyhedralDomain.h"
#include "Statement.h"

#include "core/SetRelation.h"
#include "core/Hadron.h"
#include "core/Descriptor.h"
#include "core/Domain.h"
#include "core/util/DebugStream.h"

//#include "HierarGraph.h"
#include "HierarPart.h"

#include "core/util/Time.h"

//#include <cloog/cloog.h>

namespace AlgoTrans {
	template <class R> class CComputation;
	template <class R> class CStatement;
	template <class R> class CReference;

	template <class R> class CModule;
	template <class R> class CCone;
	template <class R> class CAffineTransformation;
	template <class R> class CLattice;

	template <class R> class CVariable {
		template <class G> friend void CComputation<G>::registerVariable(CVariable<G>& variable);
		template <class G> friend void CComputation<G>::registerReference(CReference<G>& reference);
		typedef CReference<R> Reference;
	private:
		CComputation<R>& computation;
		vector<Reference*> references;
		int index;

		//CDisjunctiveDomain<Set>* dataSpace;

		CVariable();

	public:
		CVariable(CComputation<R>& iComputation);

		Reference& getReference(int index) { return *references[index]; }
		void registerReference(Reference& reference);
		int getReferenceCount() const { return references.size(); }

		int getIndex() const { return index; }
		std::string getName() const { return "v" + CInteger(getIndex()).toString(); }

		//const CDisjunctiveDomain<Set>& getDataSpace() const { return *dataSpace; }
	};

	enum UReadWrite { READ, WRITE };

	template <class R> class CReference {
		template <class G> friend void CStatement<G>::registerReference(CReference<G>& reference);
		template <class G> friend void CVariable<G>::registerReference(CReference<G>& reference);
		template <class G> friend void CComputation<G>::registerReference(CReference<G>& reference);
		typedef CComputation<R> Computation;
		typedef CStatement<R> Statement;
		typedef CVariable<R> Variable;
	private:
		Computation& computation;
		Statement& statement;
		Variable& variable;
		UReadWrite readWrite;
		CAffineTransformation<R> transformation;

		int indexInStatement;
		int indexInVariable;
	public:
		CReference(Computation& iComputation, Statement& iStatement, Variable& iVariable,
				   UReadWrite iReadWrite, const CAffineTransformation<R>& iTransformation);

		bool isWriteReference() const { return readWrite == WRITE; }
		Statement& getStatement() { return statement; };
		const Variable& getVariable() const { return variable; };
		const CAffineTransformation<R>& getTransformation() const { return transformation; };

		void print() const { printf("%s\n", toString().c_str()); };
		std::string toString() const;
		std::string toStringHtml() const;
	};


	template <class R> class CComputation {
		typedef CStatement<R> Statement;
		typedef CVariable<R> Variable;
		typedef CReference<R> Reference;

		typedef DeskriptorSet<Deskriptor<R> > CDeskriptorSet;
		typedef Hadron<DeskriptorSet<Deskriptor<R> > > CHadron;
		typedef Flat<Hadron<DeskriptorSet<Deskriptor<R> > > > Set;
		typedef SetRelation<Set> SetRel;
		typedef CDisjunctiveDomain<Set> DDomain;
		typedef SetRelation<DDomain> DomainRelation;
		typedef CGraph<Statement*, DomainRelation> RDG;
		typedef CVertex<Statement*, DomainRelation> RDGVertex;
		typedef CEdge<Statement*, DomainRelation> RDGEdge;

		friend class CStatement<R>;

		//typedef CHierarNodeLeafIterator<CComputation, CPolyhedralDomain<R>, Statement> ExecutionLeafIterator;
		typedef CHierarPart<CComputation, R> ExecutionOrderHP;
		typedef CHierarLeafPart<CComputation, R> ExecutionOrderHLP;
		typedef CHierarFinitePart<CComputation, R> ExecutionOrderHFP;
		typedef CHierarSetPart<CComputation, R> ExecutionOrderHSP;
		typedef CHierarScatterPart<CComputation, R> ExecutionOrderHCP;

		//typedef CHierarNode<CComputation, CPolyhedralDomain<R>, Statement> ExecutionOrderHP;

		//typedef typename CLatticeRelation<R>::T LatticeRelation;
		//typedef typename CRADGVertex<CStatement<R>*, LatticeRelation>::T CRADGLatticeVertex;

		//typedef CRADGraph<CStatement<R>*, SetRelation<CFlat<CLattice<R> > > > LatticeRADGraph;
		//typedef CRADGraph<CStatement<R>*, SetRelation<CFlat<CModule<R> > > > ModuleRADGraph;

	private:
		vector<Statement*> statements;
		vector<Variable*> variables;
		vector<Reference*> references;

		std::vector<std::vector<Set*> > cachedDomainCones;

		ExecutionOrderHP* originalExecutionOrder;
		bool cachedDomainsInitialized;

		//template <class SetClass> CGraph<Statement*, SetRelation<CFlat<SetClass> > > calculateApproximatedRSDG(bool polyhedral);
		CGraph<DDomain, SetRel > calculateSpace_RDG();
		RDG getRDG();
		CGraph<DDomain, SetRel > calculateTime_RDG();
		CGraph<Statement*, SetRel>* calculateDependenceClosure();

		Set getAffineRelationDescriptorsFromReferencePair(Reference& refFrom, Reference& refTo);
		Set getConvexRelationDescriptorsFromReferencePair(Reference& refFrom, Reference& refTo);
		DDomain getRelationDescriptorsFromReferencePair(Reference& refFrom, Reference& refTo);
		DDomain getReflexiveRelationDescriptorsFromReferencePair(Reference& refFrom, Reference& refTo);

		typedef CHierarPart<CComputation, R> ExecutionOrderHN;
	public:
		ExecutionOrderHN* executionOrder;
		typedef CGraph<Statement*, SetRel> AbsRDG;
		typedef CVertex<Statement*, SetRel> AbsRDGVertex;
		typedef CEdge<Statement*, SetRel> AbsRDGEdge;

		CHadron* spaceSolutions;
		CHadron* timeSolutions;

		SetRelation<Set> getSetRelationFromReferencePair(Reference& refFrom, Reference& refTo, bool polyhedral);
		SetRelation<DDomain> getDisjunctiveSetRelationFromReferencePair(Reference& refFrom, Reference& refTo);
		SetRelation<DDomain> getReflexiveDisjunctiveSetRelationFromReferencePair(Reference& refFrom, Reference& refTo);

		Set getReferenceSpace(Reference& fromRef, Reference& toRef);
		CDisjunctiveDomain<Set> getDependenceCone(Reference& fromRef, Reference& toRef);
		CDisjunctiveDomain<Set> getReflexiveDependenceCone(Reference& fromRef, Reference& toRef);
		CDisjunctiveDomain<Set> getDomainCone(Statement& firstStatement, Statement& secondStatement);
		CDisjunctiveDomain<Set> getPairDomain(Statement& firstStatement, Statement& secondStatement);
		CDisjunctiveDomain<Set> getLexiCone(Statement& firstStatement, Statement& secondStatement);
		CDisjunctiveDomain<Set> getLexiCone(int parameterCount, int commonSpaceDim, int firstDim, int secondDim, bool firstSyntPrecedesSecond);
		bool doesSyntacticallyPrecede(Statement& firstStatement, Statement& secondStatement) const;

		//ModuleRADGraph* dependenceFlatClosure;
		//LatticeRADGraph* dependenceLatticeClosure;

		DebugOutStream *osHtml;

		void initCachedDomains();

		void registerReference(Reference& reference);
		void registerVariable(Variable& variable);
		void registerStatement(Statement& statement);

		//void setOriginalExecutionOrder(CHierarNode<CComputation, R>* newOriginalExecutionOrder) { originalExecutionOrder = newOriginalExecutionOrder; }
		void setOriginalExecutionOrder(CHierarPart<CComputation, R>* newOriginalExecutionOrder) { originalExecutionOrder = newOriginalExecutionOrder; }

		int getParameterCount() const;

		CComputation() : originalExecutionOrder(NULL), cachedDomainsInitialized(false) { };

		Statement& addStatement() { return *(new CStatement<R>(*this)); }
		Variable& addVariable();
		Reference& addReference(Statement& statement, Variable& variable, UReadWrite readWrite, CAffineTransformation<R> transformation);

		int getStatementCount() const { return statements.size(); }
		int getVariableCount() const { return variables.size(); }

		const Statement& getStatement(int index) const { return *statements[index]; }
		const Variable& getVariable(int index) const;

		template <class PartitionClass>
		void parallellize();

		CHadron calcMultiTupleRelation(AbsRDG& rdg);
		bool nonPartitioningMultiTupleRelation(RDG& rdg, CHadron& mtr);
		void computeImplementation();

		CHierarPart<CComputation<R>, R>* computeConstantSpacePartitions(RDG& rdg, bool nextConstant);
		CHierarPart<CComputation<R>, R>* computeConstantTimePartitions(RDG& rdg);
		CHierarPart<CComputation<R>, R>* computeAffineSpacePartitions(RDG& rdg);
		CHierarPart<CComputation<R>, R>* computeAffineTimePartitions(RDG& rdg);
		bool noDependenciesLeft(RDG& rdg);

		void dismissDependencies(RDG* rdg, CHadron partitions);
		//CHierarSetPart<CComputation<R>, R>* scatterDomains(RDG& rdg, std::vector<Set> scatterFuncs);
		//CHierarSetPart<CComputation<R>, R>* scatterDomainsByCompoundVector(RDG& rdg, CHadron scatterFuncs);

		AbsRDG convertRDGToSpaceRDG(RDG& rdg);
		AbsRDG convertRDGToTimeRDG(RDG& rdg);

		void performAffineSetSpacePartitioning();
		//void performAffineLatticeSpacePartitioning();

		void analyseCommunicationPartitions();

		string toStringHtml();
	};

}

#include "Variable.cpp"
#include "Reference.cpp"

namespace AlgoTrans {
	/*template <class R> template <class SetClass>
	CGraph<CStatement<R>*, SetRelation<CFlat<SetClass> > >
	CComputation<R>::calculateApproximatedRSDG(bool polyhedral) {
		typedef SetRelation<CFlat<SetClass> > SetRelation;
		typedef CGraph<CStatement<R>*, SetRelation> RADGraph;
		typedef CVertex<CStatement<R>*, SetRelation> RADGVertex;

		if (osHtml != NULL) *osHtml << "Parameters:" << getParameterCount() << "<br/><b>Basic edge relations:</b><br/><table>"
			<< "<tr><td>Variable</td><td>Left statement</td><td>Right statement</td><td>Pairwise relation</td><td>Symmetric relation</td></tr>";

		RADGraph radg = RADGraph();
		for (unsigned int q = 0; q < statements.size(); q++) radg.addVertex(statements[q]);
		radg.osHtml = osHtml;

		for (unsigned int a = 0; a < variables.size(); a++) { Variable& variable = *variables[a];
			// We will want to reduce the list to unique (statement, transformation) references for performance reasons (later)
			for (int r = 0; r < variable.getReferenceCount(); r++) { const Reference& refR = variable.getReference(r);
				RADGVertex& resultVxV = radg(refR.getStatement().getIndex());

				for (int s = 0; s <= r; s++) { const Reference& refS = variable.getReference(s);
					RADGVertex& resultVxW = radg(refS.getStatement().getIndex());

					if (osHtml != NULL) *osHtml << "<tr><td>" << a << "</td><td>" << refR.getStatement().getIndex() << "</td><td>" << refS.getStatement().getIndex() << "</td>";

					SetRelation mLR = getSetRelationFromReferencePair<SetClass>(refR, refS, polyhedral);
					SetRelation mRL = getSetRelationFromReferencePair<SetClass>(refS, refR, polyhedral);

					// XXX: Project away the kernel of the iteration space (? or keep it, and just be "aware" of it?)

					radg.addEdge(resultVxV, resultVxW, mLR);
					radg.addEdge(resultVxW, resultVxV, mRL);

					if (osHtml != NULL) *osHtml << "<td>" << mLR.toStringHtml() << "</td><td>" << mRL.toStringHtml() << "</td></tr>";
				}
			}
		}

		if (osHtml != NULL) *osHtml << "</table>";

		return radg;
	}*/

	template <class R>
	bool
	CComputation<R>::noDependenciesLeft(RDG& rdg) {
		ASSERT(false); // deprecated
		return false;
		// Check whether recursion is complete
/*		if (rdg.getVertexCount() == 1) {
			// Check whether dimensionality of self-dependence == 0
			vector<RDGEdge* > edges = rdg.getEdgesBetween(rdg(0), rdg(0));
			for (unsigned int i = 0; i < edges.size(); i++) {
				DomainRelation sr = edges[i]->getData();
				sr.print();
				DeskriptorSet<R> cp = DeskriptorSet<R>::unitDeskriptorSet(getParameterCount() + 1, true)
				                   << DeskriptorSet<R>::zeroDeskriptorSet(edges[i]->getData().getSpaceDimension() - getParameterCount(), true);
				cP.print();
				DomainRelation cpsr = sr * Hadron(C, cp);
				cpsr.print();

				int q = edges[i]->getData().getSpaceDimension() - getParameterCount();
				DeskriptorSet<R> ci = DeskriptorSet<R>::zeroDeskriptorSet(getParameterCount() + 1, true)
				                   << DeskriptorSet<R>::unitDeskriptorSet(q/2, true)
				                   << -DeskriptorSet<R>::unitDeskriptorSet(q/2, true);


				int b = cpsr.getDimensionality();

				DomainRelation cpsrci = cpsr && Hadron(C, ci);
				cpsrci.print();
				int c = cpsrci.print();

				CInteger(b).print();
				CInteger(c).print();
				CInteger(rdg(0).getData().getDimensionality()).print();
				if (edges[i]->getData().getDimensionality() == rdg(0).getData().getDimensionality() > 0) return false;
			}
			return true;
		} else return false;*/
	}

	template <class R>
	CHierarPart<CComputation<R>, R>*
	CComputation<R>::computeConstantSpacePartitions(RDG& rdg, bool nextConstant) {

		H("<h2>Constant space partitioning</h2>");
		H("<b>Component ? </b><br />");
		//printf("CS\n");
		//RDG* gJ = rdg.joinConcurrentEdges();

		AbsRDG rdgSpace = convertRDGToSpaceRDG(rdg);
		rdgSpace.osHtml = osHtml;

		H(Html::startTable("Common Data Access Graph", rdg.getVertices().size()));
		CGraph<RDGVertex*, bool> levelCDAG = CGraph<RDGVertex*, bool>(); // common data access graph
		ITT(vector<RDGVertex*>, v, rdg.getVertices()) levelCDAG.addVertex(*v);
		ITT(vector<RDGVertex*>, v, rdg.getVertices()) { H("<tr>");
			ITT(vector<RDGVertex*>, w, rdg.getVertices()) {
				vector<AbsRDGEdge*> edges = rdgSpace.getEdgesBetween((**v).getIndex(), (**w).getIndex());
				bool empty = true;
				if (edges.size() > 0) {
					//edges[0]->getData().print();
					//Bool(edges[0]->getData().isEmpty()).print();
					if (!(edges[0]->getData().isEmpty())) {
						levelCDAG.addEdge(levelCDAG((*v)->getIndex()), levelCDAG((*w)->getIndex()), true);
						empty = false;
					}
				}
				//H("<td>");
				//ITT(vector<RDGEdge*>, e, edges) { H((*e)->getData().toStringHtml()); }
				//H("</td>");
				H(std::string("<td><center>") + (empty ? "o" : "x") + "</center></td>");
			} H("</tr>");
		} H("</table>");

		H(levelCDAG.saveToSVGString(svgCount++));

		H(Html::startTable("Resulting Components", 1));

		vector<vector<RDGVertex* > > connComps = levelCDAG.getConnectedComponentsData();
		vector<vector<AbsRDGVertex* > > connCompsAbs = vector<vector<AbsRDGVertex* > >();
		for (unsigned int q = 0; q < connComps.size(); q++) {
			vector<AbsRDGVertex* > vec = vector<AbsRDGVertex* >();
			for (unsigned int r = 0; r < connComps[q].size(); r++) {
				vec.push_back(&(rdgSpace.getVertex((connComps[q])[r]->getIndex())));
			}
			connCompsAbs.push_back(vec);
		}
		vector<AbsRDG> subgraphs = rdgSpace.getVertexSetInducedSubGraphs(connCompsAbs);
		vector<RDG> subrdgs = rdg.getVertexSetInducedSubGraphs(connComps);

		H("<h2>Affine space partitioning</h2>");
		vector<ExecutionOrderHP*> subHPs = vector<ExecutionOrderHP*>();
		for (unsigned int q = 0; q < subrdgs.size(); q++) {
				resumeTimer(tcComputationClosureSpace, "Computation: Space Closure");
				//CInteger(subgraphs.size()).print();
				//CInteger(subgraphs[q].getVertexCount()).print();
				AbsRDG* rdugSpaceClosed = subgraphs[q].calcJoinTransitiveClosure();
				pauseTimer(tcComputationClosureSpace);

				CHadron spaceSolutions = calcMultiTupleRelation(*rdugSpaceClosed);
				spaceSolutions.print();

				// Check whether this is the last partition
				ExecutionOrderHP* subPart;
				if (subrdgs[q].getVertexCount() == 1) {
					spaceSolutions.getDeskriptors(G).print();
					subrdgs[q](0).getData()->getFullIterationSpace().print();
					CInteger(subrdgs[q](0).getData()->getIterationSpaceDim()).print();
					if (spaceSolutions.getDeskriptors(G).getSize() == 1 + getParameterCount() + subrdgs[q](0).getData()->getIterationSpaceDim()) {
						subPart = new ExecutionOrderHLP(*rdg(0).getData());
					}
				} else if (nextConstant) {
					subPart = computeConstantTimePartitions(subrdgs[q]);
				} else {
					subPart = computeAffineTimePartitions(subrdgs[q]);
				}
				if (nonPartitioningMultiTupleRelation(subrdgs[q], spaceSolutions)) {
					subHPs.push_back(subPart);
				} else {
					subHPs.push_back(ExecutionOrderHCP::fromUnorderedSubNode(spaceSolutions, subPart));
				}
		} H("</table>");

		if (subHPs.size() > 1) { return ExecutionOrderHFP::fromUnorderedSubNodes(subHPs); }
		else { return subHPs[0]; }
	}

	template <class R>
	typename CComputation<R>::AbsRDG
	CComputation<R>::convertRDGToSpaceRDG(RDG& rdg) {
		AbsRDG result = AbsRDG(); // common data access graph

		H(Html::startTable("Converting RDG to Space RDG", rdg.getVertices().size()));
		ITT(vector<RDGVertex*>, v, rdg.getVertices()) result.addVertex((*v)->getData());
		ITT(vector<RDGVertex*>, v, rdg.getVertices()) { AbsRDGVertex& vv = result.getVertex((*v)->getIndex());
			ITT(vector<RDGVertex*>, w, rdg.getVertices()) { AbsRDGVertex& ww = result.getVertex((*w)->getIndex()); H("<tr>");
				vector<RDGEdge*> edges = rdg.getEdgesBetween(**v, **w);
				DomainRelation* relation = NULL;
				bool empty = true;
				ITT(vector<RDGEdge*>, e, edges) if (*e != NULL) {
					if (empty) { relation = new DomainRelation((*e)->getData()); }
					else { *relation = ((*relation) || (*e)->getData()); }
					empty = false;
				}
				if (relation != NULL) {
				    if (!empty) result.addEdge(vv, ww, SetRel(relation->getConstDimension(), relation->getPDimension(), relation->getQDimension(), relation->getBidirectionalHull()));
				    H(std::string("<td>") + CInteger((*v)->getIndex()).toString() + "</td><td>" + CInteger((*w)->getIndex()).toString() + "</td><td>"
				      + SetRel(relation->getConstDimension(), relation->getPDimension(), relation->getQDimension(), relation->getBidirectionalHull()).toStringHtml() + "</center></td>");
				    delete relation;
				}
				H("</tr>");
			}
		} H("</table>");

		return result;
	}

	template <class R>
	typename CComputation<R>::AbsRDG
	CComputation<R>::convertRDGToTimeRDG(RDG& rdg) {
		AbsRDG result = AbsRDG(); // common data access graph

		H(Html::startTable("Converting RDG to Time RDG", rdg.getVertices().size()));
		ITT(vector<RDGVertex*>, v, rdg.getVertices()) result.addVertex((*v)->getData());
		ITT(vector<RDGVertex*>, v, rdg.getVertices()) { AbsRDGVertex& vv = result.getVertex((*v)->getIndex()); H("<tr>");
			ITT(vector<RDGVertex*>, w, rdg.getVertices()) { AbsRDGVertex& ww = result.getVertex((*v)->getIndex());
				vector<RDGEdge*> edges = rdg.getEdgesBetween(**v, **w);
				DomainRelation* relation = NULL;
				bool empty = true;
				ITT(vector<RDGEdge*>, e, edges) if (*e != NULL) {
					if (empty) { relation = new DomainRelation((*e)->getData()); }
					else { *relation = ((*relation) || (*e)->getData()); }
					empty = false;
				}
				if (relation != NULL) {
				    if (!empty) result.addEdge(vv, ww, SetRel(relation->getConstDimension(), relation->getPDimension(), relation->getQDimension(), relation->getUnidirectionalHull()));
				    H(std::string("<td><center>") + SetRel(relation->getConstDimension(), relation->getPDimension(), relation->getQDimension(), relation->getUnidirectionalHull()).toString() + "</center></td>");
				    delete relation;
				}
			} H("</tr>");
		} H("</table>");

		return result;
	}

	/*template <class R>
	CHierarSetPart<CComputation<R>, R>*
	CComputation<R>::scatterDomains(RDG& rdg, std::vector<Set> scatterFuncs) {
        DDomain newDom = DDomain(scatterFuncs.size());
		ITT(vector<RDGVertex*>, v, rdg.getVertices()) {
			DDomain& d = rdg.getData();

			for (int q = 0; q < d.getSize(); q++) {

			}
		}
	}

	template <class R>
	CHierarSetPart<CComputation<R>, R>*
	CComputation<R>::scatterDomainsByCompoundVector(RDG& rdg, CHadron scatterFuncs) {
		std::vector<Set> scatterFuncsV;
		int offset = 0;
		ITT(vector<RDGVertex*>, v, rdg.getVertices()) {
			DDomain& d = rdg.getData();

			std::vector<int> retainedDimensions;
			for (int q = 0; q < d.getSpaceDimension(); q++) retainedDimensions.push_back(offset + q);
			Set scatterF = scatterFuncs.getProjection(retainedDimensions);
			scatterF.print();
            scatterFuncsV.push_back(scatterF);

            offset += d.getSpaceDimension();
		}
	}*/

	template <class R>
	void CComputation<R>::dismissDependencies(RDG* rdg, CHadron partitions) {
		partitions.print();

		DeskriptorSet<Deskriptor<R> > g = partitions.getDeskriptors(G);

		ITT(RDG, v, *rdg) {
			ITT(RDG, w, *rdg) {
				int constIstart = 0;
				int constJstart = 0;
				int dimIstart = 0;
				int dimJstart = 0;
				int offset = 0;
				//CInteger(rdg.getVertexCount()).print();
				for (int q  = 0; q < rdg->getVertexCount(); q++) {
					if (q == (*v)->getIndex()) {
						constIstart = offset;
						dimIstart = 1 + getParameterCount() + constIstart;
					}
					if (q == (*w)->getIndex()) {
						constJstart = offset;
						dimJstart = 1 + getParameterCount() + constJstart;
					}
					offset += 1 + rdg->getVertex(q).getData()->getFullIterationSpace().getSpaceDimension();
				}

				vector<RDGEdge*> edges = rdg->getEdgesBetween((*v)->getIndex(), (*w)->getIndex());

				vector<int> retainedDimensions = vector<int>();
				for (int q = 0; q < (*v)->getData()->getFullIterationSpace().getSpaceDimension() + 1; q++) {
					retainedDimensions.push_back(constIstart + q);
				}
				if ((*v)->getIndex() != (*w)->getIndex()) {
					for (int q = 0; q < (*w)->getData()->getFullIterationSpace().getSpaceDimension() + 1; q++) {
						retainedDimensions.push_back(constJstart + q);
					}
				}
				DeskriptorSet<Deskriptor<R> > pairPartition = g.getGeneratingProjection(retainedDimensions);
				if ((*v)->getIndex() == (*w)->getIndex()) pairPartition = pairPartition << pairPartition;
				pairPartition.print();

				DeskriptorSet<Deskriptor<R> > constDeps = DeskriptorSet<Deskriptor<R> >(1 + getParameterCount() + (*v)->getData()->getIterationSpaceDim() + (*w)->getData()->getIterationSpaceDim());
				for (int q = 0; q < pairPartition.getSize(); q++) {
					typedef typename Deskriptor<R>::Vector Vector;
					Vector vec = Vector();
					pairPartition[q].print();
					for (int i = 0; i < getParameterCount() + 1; i++) {
						vec.appendElement(pairPartition[q][constIstart + i] - pairPartition[q][constJstart + i]);
					}
					for (int i = 0; i < (*v)->getData()->getIterationSpaceDim(); i++) {
						vec.appendElement(pairPartition[q][dimIstart + i]);
					}
					for (int j = 0; j < (*w)->getData()->getIterationSpaceDim(); j++) {
						vec.appendElement(-pairPartition[q][dimJstart + j]);
					}
					vec.print();
					if (!vec.isZero()) constDeps.addDeskriptor(Deskriptor<R>(vec, true));
				}

				ITT(vector<RDGEdge*>, e, edges) {
					DomainRelation r = (*e)->getData();
					r.print();
					for (int s = 0; s < r.getElementCount(); s++) {
						r[s].print();
						r[s] = r[s] * CHadron(C, constDeps);
						r[s].print();
					}
					r.print();
				}
			}
		}

	}

	template <class R>
	CHierarPart<CComputation<R>, R>*
	CComputation<R>::computeAffineSpacePartitions(RDG& rdg) {
		H("<h2>Affine space partitioning</h2>");
		typedef CHierarPart<CComputation, R> ExecutionOrderHN;

		printf("AS\n");

		AbsRDG rdgSpace = convertRDGToSpaceRDG(rdg);
		rdgSpace.osHtml = osHtml;
		//rdgSpace.print();

		typedef typename SetRel::DisjunctiveHull_CombinationOperation CombOp;
		typedef typename SetRel::template Operations<CombOp> Ops;
		Ops ops = Ops();
		//RDG* rdugSpace = rdgSpace.combineConcurrentEdges(ops);
		//rdgSpace.getUniqueEdgeBetween(rdgSpace.getVertex(0), rdgSpace.getVertex(0)).getData().print();
		//rdugSpace->osHtml = osHtml;

		resumeTimer(tcComputationClosureSpace, "Computation: Space Closure");
		AbsRDG* rdugSpaceClosed = rdgSpace.calcJoinTransitiveClosure();
		pauseTimer(tcComputationClosureSpace);

		CHadron spaceSolutions = calcMultiTupleRelation(*rdugSpaceClosed);
		spaceSolutions.print();

		//scatterDomainsByCompoundVector(rdg, *spaceSolutions);
		// Update rdg

		//dismissDependencies(&rdg, spaceSolutions);

		return ExecutionOrderHCP::fromUnorderedSubNode(spaceSolutions, computeConstantTimePartitions(rdg));
		/*typedef SetRelation<Set> SetRelation;

		H("<b>Component " << "?" << "</b><br />");

		CGraph<CStatement<R>*, SetRelation<Set> >  gJ = g.joinConcurrentEdges();

		RDG rdgSpace = calculateSpace_RDG();
		rdgSpace.osHtml = osHtml;

		g.calcJoinTransitiveClosure();*/

		// Extract partitions

		// Generate LevelData: PolyhedralDomain for resulting hp

		// Reduce relation of g to reflect extracted partitions
		//if (osHtml != NULL) levelCDAG.saveToSVGStream(*osHtml);

		//H(Html::startTable("Resulting Components", 1));
		//vector<RADGraph> subgraphs = g.getVertexSetInducedSubGraphs(levelCDAG.getConnectedComponents());
		/*ITT(EdgeIt, e, subGraph) {
			RDVertex& from = (*e)->getFromVertex();
			RDVertex& to = (*e)->getToVertex();
			// XXX : Reduce the corresponding relation with the chosen partitions
			(*e)->setData((*e)->getData());
		} H("</table>");*//*

		//return ExecutionOrderHP::fromUnorderedSubNodes(pd, computeConstantSpacePartitions(subGraph));*/
	}

	template <class R>
	CHierarPart<CComputation<R>, R>*
	CComputation<R>::computeAffineTimePartitions(RDG& rdg) {
		//if (noDependenciesLeft(rdg)) { return new ExecutionOrderHLP(*rdg(0).getData()); }

		H("<h2>Affine time partitioning</h2>");
		typedef CHierarPart<CComputation, R> ExecutionOrderHN;

		AbsRDG rdgSpace = convertRDGToTimeRDG(rdg);
		rdgSpace.osHtml = osHtml;
		//rdgSpace.print();

		typedef typename SetRel::DisjunctiveHull_CombinationOperation CombOp;
		typedef typename SetRel::template Operations<CombOp> Ops;
		Ops ops = Ops();
		//RDG* rdugSpace = rdgSpace.combineConcurrentEdges(ops);
		//rdgSpace.getUniqueEdgeBetween(rdgSpace.getVertex(0), rdgSpace.getVertex(0)).getData().print();
		//rdugSpace->osHtml = osHtml;

		resumeTimer(tcComputationClosureSpace, "Computation: Time Closure");
		AbsRDG* rdugSpaceClosed = rdgSpace.calcJoinTransitiveClosure();
		pauseTimer(tcComputationClosureSpace);

		CHadron spaceSolutions = calcMultiTupleRelation(*rdugSpaceClosed);
		spaceSolutions.print();

		//scatterDomainsByCompoundVector(rdg, *spaceSolutions);
		// Update rdg

		dismissDependencies(&rdg, spaceSolutions);

		return ExecutionOrderHCP::fromLexicoOrderedSubNode(spaceSolutions, computeConstantSpacePartitions(rdg, true));
		/*typedef SetRelation<Set> SetRelation;

		H("<b>Component " << "?" << "</b><br />");

		CGraph<CStatement<R>*, SetRelation<Set> >  gJ = g.joinConcurrentEdges();

		RDG rdgSpace = calculateSpace_RDG();
		rdgSpace.osHtml = osHtml;

		g.calcJoinTransitiveClosure();*/

		// Extract partitions

		// Generate LevelData: PolyhedralDomain for resulting hp

		// Reduce relation of g to reflect extracted partitions
		//if (osHtml != NULL) levelCDAG.saveToSVGStream(*osHtml);

		//H(Html::startTable("Resulting Components", 1));
		//vector<RADGraph> subgraphs = g.getVertexSetInducedSubGraphs(levelCDAG.getConnectedComponents());
		/*ITT(EdgeIt, e, subGraph) {
			RDVertex& from = (*e)->getFromVertex();
			RDVertex& to = (*e)->getToVertex();
			// XXX : Reduce the corresponding relation with the chosen partitions
			(*e)->setData((*e)->getData());
		} H("</table>");*//*

		//return ExecutionOrderHP::fromUnorderedSubNodes(pd, computeConstantSpacePartitions(subGraph));*/
	}

	template <class R>
	CHierarPart<CComputation<R>, R>*
	CComputation<R>::computeConstantTimePartitions(RDG& rdg) {
		//if (noDependenciesLeft(rdg)) { return new ExecutionOrderHLP(*rdg(0).getData()); }

		H("<h2>Constant time partitioning</h2>");
		typedef vector<CVertex<RDG, bool>* > GraphVertexVector;
		typedef CIxLink<RDGVertex, CVertex<ExecutionOrderHP*, Bool>*, vector<CVertex<ExecutionOrderHP*, Bool>* > > VertIxLink;

		H("<b>Component ? </b><br />");

		//RDG* gJ = rdg.joinConcurrentEdges();

		H(Html::startTable("Common Data Access Graph", rdg.getVertices().size()));
		CGraph<RDGVertex*, bool> levelCDAG = CGraph<RDGVertex*, bool>(); // common data access graph
		ITT(vector<RDGVertex*>, v, rdg.getVertices()) levelCDAG.addVertex(*v);
		ITT(vector<RDGVertex*>, v, rdg.getVertices()) { H("<tr>");
			ITT(vector<RDGVertex*>, w, rdg.getVertices()) {
				vector<RDGEdge*> edges = rdg.getEdgesBetween(**v, **w);
				//printf("c\n");
				bool empty = true;
				ITT(vector<RDGEdge*>, e, edges) {
					if (*e != NULL) (*e)->getData().print();
					if ((*e != NULL) && (!(*e)->getData().isEmpty())) { empty = false; break; }
				}
					H("<td>");
					ITT(vector<RDGEdge*>, e, edges) { H((*e)->getData().toStringHtml()); }
					H("</td>");
				  if (!empty) levelCDAG.addEdge(levelCDAG((*v)->getIndex()), levelCDAG((*w)->getIndex()), true);
				  H(std::string("<td><center>") + (empty ? "o" : "x") + "</center></td>");
			} H("</tr>");
		} H("</table>");

		H(levelCDAG.saveToSVGString(svgCount++));

		H(Html::startTable("Resulting Components", 1));

		CGraph<RDG, bool> sccGraph = rdg.getVertexSetGraphInducedGraph(levelCDAG.getStronglyConnectedDecompositionSCCsData());
		vector<ExecutionOrderHP*> subHPs = vector<ExecutionOrderHP*>();
		ITT(GraphVertexVector, hI, sccGraph) {
			ITT(vector<RDGVertex*>, v, ((*hI)->getData()).getVertices()) { if (v != (*hI)->getData().getVertices().begin()) H(", ");
				H(CInteger((*v)->getIndex()).toString());
				//ExecutionOrderHN resultHN = ExecutionOrderHN(*this, CDisjunctiveDomain<Set>::universe(parentHN->getLevelData()->getSpaceDimension()));
			}
			subHPs.push_back(computeConstantSpacePartitions((*hI)->getData(), false));
		} H("</table>");

		if (subHPs.size() == 1) return subHPs[0];

		ExecutionOrderHFP* hfp = ExecutionOrderHFP::fromUnorderedSubNodes(subHPs);

		VertIxLink vLink = VertIxLink(hfp->cellGraph.getVertices());
		typedef CVertex<RDG, bool> SCCVertex;
		vector<SCCVertex*> vertices = sccGraph.getVertices();
		ITT(vector<SCCVertex*>, sccA, vertices) {
			ITT(vector<SCCVertex*>, sccB, vertices) {
				typedef vector<CEdge<RDG, bool>* > EdgeVector;
				EdgeVector edges = sccGraph.getEdgesBetween((*sccA)->getIndex(), (*sccB)->getIndex());
				ITT(EdgeVector, e, edges) hfp->cellGraph.addEdge(hfp->cellGraph((*sccA)->getIndex()), hfp->cellGraph((*sccB)->getIndex()), (*e)->getData());
			}
		}

		return hfp;
	}

/*	template <class R>
	CHierarPart<CComputation<R>, R>*
	CComputation<R>::computeConstantTimePartitions(CGraph<Set, SetRelation<Set> >& g) {
		typedef SetRelation<Set> SR;
		typedef CGraph<Set, SR> RDG;
		typedef CVertex<Set, SR> RDGVertex;
		typedef CEdge<Set, SR> RDGEdge;
		//typedef typename CGraph<ExecutionOrderHP*, Bool>::iterator HPGraphVertexIterator;
		typedef vector<CVertex<RDG, bool>* > GraphVertexVector;
		typedef CIxLink<RDGVertex, CVertex<ExecutionOrderHP*, Bool>*, vector<CVertex<ExecutionOrderHP*, Bool>* > > VertIxLink;

		H("<b>Component " << "?" << "</b><br />");
		RDG* gJ = g.joinConcurrentEdges();

		H(Html::startTable("Common Data Access Graph", g.getVertices().size()));
		CGraph<RDGVertex*, bool> levelCDAG = CGraph<RDGVertex*, bool>(); // common data access graph
		ITT(vector<RDGVertex*>, v, gJ->getVertices()) levelCDAG.addVertex(*v);
		ITT(vector<RDGVertex*>, v, gJ->getVertices()) { H("<tr>");
			ITT(vector<RDGVertex*>, w, gJ->getVertices()) {
				vector<RDGEdge*> edges = gJ->getEdgesBetween(**v, **w);
				bool empty = true;
				ITT(vector<RDGEdge*>, e, edges) if ((*e != NULL) && (!(*e)->getData().isEmpty())) { empty = false; break; }
				if (!empty) levelCDAG.addEdge(levelCDAG((*v)->getIndex()), levelCDAG((*w)->getIndex()), true);
				H(std::string("<td><center>") + (empty ? "o" : "x") + "</center></td>");
			} H("</tr>");
		} H("</table>");

		if (osHtml != NULL) levelCDAG.saveToSVGStream(*osHtml);

		H(Html::startTable("Resulting Components", 1));
		CGraph<RDG, bool> sccGraph = g.getVertexSetGraphInducedGraph(levelCDAG.getStronglyConnectedDecompositionSCCs());
		vector<ExecutionOrderHP*> subHPs = vector<ExecutionOrderHP*>();
		ITT(GraphVertexVector, hI, sccGraph) {
			ITT(vector<RDGVertex*>, v, (*hI).getVertices()) { if (v != (*hI).getVertices().begin()) H(", ");
				H(CInteger((*v)->getIndex()).toString());
				//ExecutionOrderHN resultHN = ExecutionOrderHN(*this, CDisjunctiveDomain<Set>::universe(parentHN->getLevelData()->getSpaceDimension()));
			}
			subHPs.push_back(computeAffineTimePartitions(*hI));
		} H("</table>");

		ExecutionOrderHP hfp = ExecutionOrderHP::fromUnorderedSubNodes(subHPs);

		VertIxLink vLink = VertIxLink(hfp.subGraph.getVertices());
		ITT(vector<RDGVertex*>, sccA, sccGraph.getVertices()) {
			ITT(vector<RDGVertex*>, sccB, sccGraph.getVertices()) {
				typedef vector<CEdge<RDG, Bool>* > EdgeVector;
				EdgeVector edges = sccGraph.getEdgesBetween(*sccA, *sccB);
				ITT(EdgeVector, e, edges) hfp.subGraph.addEdge(vLink[*sccA], vLink[*sccB], Bool((*e)->getData()));
			}
		}

		return hfp;
	}*/

	template <class R>
	std::string CReference<R>::toString() const {
		return statement.getName() + ":" + variable.getName() + "[" + transformation.toString() + "]";
	}

	template <class R>
	std::string CReference<R>::toStringHtml() const {
		return "<table><tr><td>" + statement.getName() + "</td><td>" + variable.getName() + "</td><td>" + transformation.toStringHtml() + "</td></tr></table>";
	}


	template <class R>
	CGraph<CDisjunctiveDomain<typename CComputation<R>::Set>, SetRelation<typename CComputation<R>::Set> >
	CComputation<R>::calculateSpace_RDG() {
		if (tcComputationCalculateSpaceRDG == NULL) { tcComputationCalculateSpaceRDG = new TimeCollector("Computation: Calculate Space RDG"); timeCollectors.push_back(tcComputationCalculateSpaceRDG); }
		tcComputationCalculateSpaceRDG->resume();

		typedef SetRelation<Set> SetRelation;
		typedef CGraph<CDisjunctiveDomain<Set>, SetRelation> RDG;
		typedef CVertex<CDisjunctiveDomain<Set>, SetRelation> RDGVertex;

		if (osHtml != NULL) *osHtml << "Parameters:" << getParameterCount() << "<br/><h2>Basic edge relations:</h2><br/><table>"
			<< "<tr><td>Variable</td><td>Left statement</td><td>Right statement</td><td>Pairwise relation</td><td>Symmetric relation</td></tr>";

		RDG rdg = RDG();

		statements[0]->getFullIterationSpace();
		ITT(vector<Statement*>, s, statements) rdg.addVertex((*s)->getFullIterationSpace());
		rdg.osHtml = osHtml;

		ITT(vector<Variable*>, var, variables) {
			// PERF:  We will want to reduce the list to unique (statement, transformation) references for performance reasons (later)
			for (int r = 0; r < (*var)->getReferenceCount(); r++) { Reference& refR = (*var)->getReference(r);
				RDGVertex& resultVxV = rdg(refR.getStatement().getIndex());

				for (int s = 0; s <= r; s++) { Reference& refS = (*var)->getReference(s);
					RDGVertex& resultVxW = rdg(refS.getStatement().getIndex());

					if (osHtml != NULL) *osHtml << "<tr><td>" << (*var)->getIndex() << "</td><td>" << refR.getStatement().getIndex() << "</td><td>" << refS.getStatement().getIndex() << "</td>";
					if (osHtml != NULL) *osHtml << "<td>" << getReferenceSpace(refR, refS).toStringHtml() << "</td><td>"
							<< getDomainCone(refR.getStatement(), refS.getStatement()).toStringHtml() << "</td><td>"
							<< getLexiCone(refR.getStatement(), refS.getStatement()).toStringHtml() << "</td><td>"
							<< getReflexiveDependenceCone(refR, refS).toStringHtml() << "</td>";

					SetRelation mLR = getSetRelationFromReferencePair(refR, refS, false);
					// PERF: Just mirror the above relation (later)
					SetRelation mRL = getSetRelationFromReferencePair(refS, refR, false);

					// PERF: Project away the useless dimensions of the iteration space (? or keep them, and just be "aware" of them?)

					rdg.addEdge(resultVxV, resultVxW, mLR);
					rdg.addEdge(resultVxW, resultVxV, mRL);

					if (osHtml != NULL) *osHtml << "<td>" << mLR.toStringHtml() << "</td><td>" << mRL.toStringHtml() << "</td></tr>";
				}
			}
		}

		if (osHtml != NULL) *osHtml << "</table>";

		tcComputationCalculateSpaceRDG->pause();

		return rdg;
	}

	template <class R>
	typename CComputation<R>::RDG
	CComputation<R>::getRDG() {
		resumeTimer(tcRDG, "RDG");


		H("<h1>Retrieving RDG:</h1>"); H("Parameters:"); H(CInteger(getParameterCount()).toString());
		H("<br/><b>Basic edge relations:</b><br/><table>");
		H("<tr><td>Var</td><td>Left S</td><td>RW</td><td>Right S</td><td>RW</td><td>Ref Space</td><td>Dom Cone</td><td>Lexi Cone</td><td>Disj Set Rel</td><td>mLR</td><td>mRL</td></tr>");

		RDG rdg = RDG();

		ITT(vector<Statement*>, s, statements) rdg.addVertex(*s);
		rdg.osHtml = osHtml;

		ITT(vector<Variable*>, var, variables) {
			// PERF:  We will want to reduce the list to unique (statement, transformation) references for performance reasons (later)
			for (int r = 0; r < (*var)->getReferenceCount(); r++) { Reference& refR = (*var)->getReference(r);
				RDGVertex& resultVxV = rdg(refR.getStatement().getIndex());

				for (int s = 0; s < (*var)->getReferenceCount(); s++) { Reference& refS = (*var)->getReference(s);
					RDGVertex& resultVxW = rdg(refS.getStatement().getIndex());

					H("<tr><td>"); H(CInteger((*var)->getIndex()).toString());
					H("</td><td>"); H(CInteger(refR.getStatement().getIndex()).toString());
					H("</td><td>"); H(refR.isWriteReference() ? "W" : "R");
					H("</td><td>"); H(CInteger(refS.getStatement().getIndex()).toString()); H("</td>");
					H("</td><td>"); H(refS.isWriteReference() ? "W" : "R");
					H("<td>"); H(getReferenceSpace(refR, refS).toStringHtml()); H("</td><td>");
							H(getDomainCone(refR.getStatement(), refS.getStatement()).toStringHtml()); H("</td><td>");
							H(getLexiCone(refR.getStatement(), refS.getStatement()).toStringHtml()); H("</td><td>");
							H(getDisjunctiveSetRelationFromReferencePair(refR, refS).toStringHtml()); H("</td>");

					//printf("%i -- %i\n", refR.getStatement().getIndex(), refS.getStatement().getIndex());
					//getLexiCone(refR.getStatement(), refS.getStatement()).print();
					DomainRelation mLR = getDisjunctiveSetRelationFromReferencePair(refR, refS);
					//mLR.print();
					// PERF: Just mirror the above relation (later)
					DomainRelation mRL = getDisjunctiveSetRelationFromReferencePair(refS, refR);
					//mRL.print();

					// PERF: Project away the useless dimensions of the iteration space (? or keep them, and just be "aware" of them?)

					rdg.addEdge(resultVxV, resultVxW, mLR);
					rdg.addEdge(resultVxW, resultVxV, mRL);

					H("<td>"); H(mLR.toStringHtml()); H("</td><td>"); H(mRL.toStringHtml()); H("</td></tr>");
				}
			}
		} H("</table>");

		pauseTimer(tcRDG);

		return rdg;
	}


	template <class R>
	CGraph<CDisjunctiveDomain<typename CComputation<R>::Set>, SetRelation<typename CComputation<R>::Set> >
	CComputation<R>::calculateTime_RDG() {
		typedef SetRelation<Set> SetRelation;
		typedef CGraph<CDisjunctiveDomain<Set>, SetRelation> RDG;
		typedef CVertex<CDisjunctiveDomain<Set>, SetRelation> RDGVertex;

		if (tcComputationCalculateTimeRDG == NULL) { tcComputationCalculateTimeRDG = new TimeCollector("Computation: Calculate Time RDG"); timeCollectors.push_back(tcComputationCalculateTimeRDG); }
		tcComputationCalculateTimeRDG->resume();

		if (osHtml != NULL) *osHtml << "Parameters:" << getParameterCount() << "<br/><b>Basic edge relations:</b><br/><table>"
			<< "<tr><td>Var</td><td>L S</td><td>R S</td><td>Left reference</td><td>Right reference</td>"
			<< "<td>Reference space</td><td>Domain Cone</td><td>LexiCone</td><td>Reflexive Dependence Cone</td>"
			<< "<td>Reflexive Dependence Cone && LexiCone</td><td>Empty</td><td>Symm Empty</td><td>Pairwise relation</td><td>Symmetric relation</td></tr>";

		RDG rdg = RDG();

		ITT(vector<Statement*>, s, statements) rdg.addVertex((*s)->getFullIterationSpace());
		rdg.osHtml = osHtml;

		ITT(vector<Variable*>, var, variables) {
			// We will want to reduce the list to unique (statement, transformation) references for performance reasons (later)
			for (int r = 0; r < (*var)->getReferenceCount(); r++) { Reference& refR = (*var)->getReference(r);
				RDGVertex& resultVxV = rdg(refR.getStatement().getIndex());

				for (int s = 0; s < (*var)->getReferenceCount(); s++) { Reference& refS = (*var)->getReference(s);
				    if (refR.isWriteReference() || refS.isWriteReference()) {
					RDGVertex& resultVxW = rdg(refS.getStatement().getIndex());

					if (osHtml != NULL) *osHtml << "<tr><td>" << (*var)->getIndex() << "</td><td>" << refR.getStatement().getIndex() << "</td><td>" << refS.getStatement().getIndex() << "</td>";
					if (osHtml != NULL) *osHtml << "<td>" << refR.toStringHtml() << "</td><td>" << refS.toStringHtml() << "</td>";
					if (osHtml != NULL) *osHtml << "<td>" << getReferenceSpace(refR, refS).toStringHtml() << "</td><td>"
							<< getDomainCone(refR.getStatement(), refS.getStatement()).toStringHtml() << "</td><td>"
							<< getLexiCone(refR.getStatement(), refS.getStatement()).toStringHtml() << "</td><td>"
							<< getReflexiveDependenceCone(refR, refS).toStringHtml() << "</td><td>"
							<< (getReflexiveDependenceCone(refR, refS) && getLexiCone(refR.getStatement(), refS.getStatement())).toStringHtml() << "</td>";

					SetRelation mLR = getSetRelationFromReferencePair(refR, refS, true);
					SetRelation mRL = getSetRelationFromReferencePair(refS, refR, true);

					if (osHtml != NULL) *osHtml << "<td>" << (((Set) mLR).isEmpty() ? "X" : "") << "</td>";
					if (osHtml != NULL) *osHtml << "<td>" << (((Set) mRL).isEmpty() ? "X" : "") << "</td>";

					// XXX: Project away the useless dimensions of the iteration space (? or keep them, and just be "aware" of them?)

					if (!(((Set) mLR).isEmpty())) rdg.addEdge(resultVxV, resultVxW, mLR);
					if (!(((Set) mRL).isEmpty())) rdg.addEdge(resultVxW, resultVxV, mRL);

					if (osHtml != NULL) *osHtml << "<td>" << mLR.toStringHtml() << "</td><td>" << mRL.toStringHtml() << "</td></tr>";
				    }
				}
			}
		}

		/*for (int v = 0; v < rdg.getVertexCount(); v++) { RDGVertex& V = rdg(v);
		    int q = V.getStatement().getIterationSpaceDim();
		    DeskriptorSet ds = DeskriptorSet::zeroDeskriptorSet(q, 1 + getParameterCount(), true);
		    SetRelation rel = SetRelation<Set>(getParameterCount(), V.getStatement().getIterationSpaceDim(), V.getStatement().getIterationSpaceDim(),
		    		              Set(1, Hadron(C, ds)));
		    rdg.addEdge(V, V, rel);
		}*/

		if (osHtml != NULL) *osHtml << "</table>";

		tcComputationCalculateTimeRDG->pause();
		return rdg;
	}

	template <class R> void CComputation<R>::registerReference(CReference<R>& reference) {
		references.push_back(&reference);
		reference.variable.registerReference(reference);
		reference.statement.registerReference(reference);
	}

	template <class R> void CComputation<R>::registerVariable(CVariable<R>& variable) {
		variable.index = variables.size();

		variables.push_back(&variable);
	}

	template <class R> void CComputation<R>::registerStatement(CStatement<R>& statement) {
		statement.index = statements.size();

		statements.push_back(&statement);
	}

	template <class R> CReference<R>& CComputation<R>::addReference(Statement& statement, Variable& variable, UReadWrite readWrite, CAffineTransformation<R> transformation) {
		return *(new CReference<R>(*this, statement, variable, readWrite, transformation));
	}


	template <class R>
	typename CComputation<R>::CHadron
	CComputation<R>::calcMultiTupleRelation(AbsRDG& rdg) {
		if (tcComputationCalcMultiTupleRelation == NULL) { tcComputationCalcMultiTupleRelation = new TimeCollector("Computation: Calculate Multi-Tuple Relation"); timeCollectors.push_back(tcComputationCalcMultiTupleRelation); }
		tcComputationCalcMultiTupleRelation->resume();

		int totalDim = 0;
		for (int i = 0; i < rdg.getVertexCount(); i++) {
			totalDim += rdg.getVertex(i).getData()->getIterationSpaceDim() + 1;
		}
		totalDim += rdg.getVertexCount() * getParameterCount();

		//Deskriptor<R> (*ZV)(int dim, bool bidir) = &(Deskriptor<R>::getZeroDeskriptor);
		Deskriptor<R> (*UV)(int dim, int unitDim, bool bidir) = &(Deskriptor<R>::getUnitDeskriptor);

		H("<b>Computing multi-tuple relation:</b>");
		H("<table><tr><td>LS</td><td>RS</td><td>rel</td><td>C,ds</td><td>hds</td><td>p</td><td>result</td></tr>");
		CHadron result = CHadron::universe(totalDim);
		for (int i = 0; i < rdg.getVertexCount(); i++) {
			for (int j = 0; j < rdg.getVertexCount(); j++) {
				const vector<AbsRDGEdge*>& edges = rdg.getEdgesBetween(i, j);
				if (edges.size() > 0) {
					CHadron rel = CHadron::source(totalDim) ^ rdg.getUniqueEdgeBetween(rdg.getVertex(i), rdg.getVertex(j)).getData();
					rel = rel.getDual();
					H("<tr><td>"); H(CInteger(i).toString()); H("</td><td>"); H(CInteger(j).toString());
					H("</td><td>"); H(rel.toStringHtml()); H("</td>");

					int constIstart = 0;
					int constJstart = 0;
					int dimIstart = 0;
					int dimJstart = 0;
					int offset = 0;
					//CInteger(rdg.getVertexCount()).print();
					for (int q  = 0; q < rdg.getVertexCount(); q++) {
						if (q == i) {
							constIstart = offset;
							dimIstart = 1 + getParameterCount() + constIstart;
						}
						if (q == j) {
							constJstart = offset;
							dimJstart = 1 + getParameterCount() + constJstart;
						}
						offset += 1 + rdg.getVertex(q).getData()->getIterationSpaceDim() + getParameterCount();
					}

					CDeskriptorSet ds = CDeskriptorSet(rel.getDimension());
					for (int q = 0; q < 1 + getParameterCount(); q++) {
						Deskriptor<R> d = -UV(rel.getDimension(), constIstart + q, true)
							        		  + UV(rel.getDimension(), constJstart + q, true)
							        		  - UV(rel.getDimension(), totalDim + q, true);
						ds.addDeskriptor(d);
					}
					for (int q = 0; q < rdg.getVertex(i).getData()->getIterationSpaceDim(); q++) {
						Deskriptor<R> d = UV(rel.getDimension(), dimIstart + q, true)
							        		  + UV(rel.getDimension(), totalDim + 1 + getParameterCount() + q, true);
						ds.addDeskriptor(d);
					}
					for (int q = 0; q < rdg.getVertex(j).getData()->getIterationSpaceDim(); q++) {
						Deskriptor<R> d = UV(rel.getDimension(), dimJstart + q, true)
							        		  - UV(rel.getDimension(), totalDim + 1 + q + rdg.getVertex(i).getData()->getFullIterationSpace().getSpaceDimension(), true);
						ds.addDeskriptor(d);
					}
					CHadron hds = CHadron(C, ds) * rel;
					H("<td>"); H(CHadron(C, ds).toStringHtml()); H("</td><td>");  H(hds.toStringHtml()); H("</td>");

					std::vector<int> retainedDims;
					for (int q = 0; q < totalDim; q++) retainedDims.push_back(q);
					CHadron p = hds.getProjection(retainedDims);
					H("<td>"); H(p.toStringHtml()); H("</td>");
					result = result * p;
					H("<td>"); H(result.toStringHtml()); H("</td></tr>");
				}
			}
	    }

		if (osHtml != NULL) *osHtml << "</table>";

		tcComputationCalcMultiTupleRelation->pause();

		return result;
	}

	template <class R>
	bool
	CComputation<R>::nonPartitioningMultiTupleRelation(RDG& rdg, CHadron& mtr) {
		if (tcComputationCalcMultiTupleRelation == NULL) { tcComputationCalcMultiTupleRelation = new TimeCollector("Computation: Calculate Multi-Tuple Relation"); timeCollectors.push_back(tcComputationCalcMultiTupleRelation); }
		tcComputationCalcMultiTupleRelation->resume();

		int totalDim = 0;

		//Deskriptor<R> (*ZV)(int dim, bool bidir) = &(Deskriptor<R>::getZeroDeskriptor);
		//Deskriptor<R> (*UV)(int dim, int unitDim, bool bidir) = &(Deskriptor<R>::getUnitDeskriptor);

		DeskriptorSet<Deskriptor<R> > g = mtr.getDeskriptors(G);
		g.print();
		for (int i = 0; i < rdg.getVertexCount(); i++) {
			for (int d = 0; d < g.getSize(); d++) { g[d].print();
				for (int j = 0; j < rdg.getVertex(i).getData()->getIterationSpaceDim(); j++) {
					if (g[d][totalDim + 1 + getParameterCount() + j] != (R) 0) return false;
				}
			}
			totalDim += rdg.getVertex(i).getData()->getIterationSpaceDim() + getParameterCount() + 1;
		}

		return true;
	}

	template <class R>
	void CComputation<R>::computeImplementation() {
		typedef CGraph<Set, SetRel> RDGSingleConjuct;

		initCachedDomains();

		RDG rdg = getRDG();
		rdg.osHtml = osHtml;

/*		AbsRDG rdgSpace = convertRDGtoSpaceRDG(rdg);
		rdgSpace.osHtml = osHtml;
		//rdgSpace.print();

		typedef typename SetRel::DisjunctiveHull_CombinationOperation CombOp;
		typedef typename SetRel::template Operations<CombOp> Ops;
		Ops ops = Ops();
		//RDG* rdugSpace = rdgSpace.combineConcurrentEdges(ops);
		//rdgSpace.getUniqueEdgeBetween(rdgSpace.getVertex(0), rdgSpace.getVertex(0)).getData().print();
		//rdugSpace->osHtml = osHtml;

		resumeTimer(tcComputationClosureSpace, "Computation: Space Closure");

		RDG* rdugSpaceClosed = rdgSpace.calcJoinTransitiveClosure();

		pauseTimer(tcComputationClosureSpace);

		spaceSolutions = new CHadron(calcMultiTupleRelation(*rdugSpaceClosed));
		spaceSolutions.print();*/

		executionOrder = computeConstantSpacePartitions(rdg, true);
#ifdef RDGDeconstruct
		RDG rdgSpace = convertRDGtoSpaceRDG(rdg);
		rdgSpace.osHtml = osHtml;
		//rdgSpace.print();

		typedef typename SetRel::DisjunctiveHull_CombinationOperation CombOp;
		typedef typename SetRel::template Operations<CombOp> Ops;
		Ops ops = Ops();
		//RDG* rdugSpace = rdgSpace.combineConcurrentEdges(ops);
		//rdgSpace.getUniqueEdgeBetween(rdgSpace.getVertex(0), rdgSpace.getVertex(0)).getData().print();
		//rdugSpace->osHtml = osHtml;

		resumeTimer(tcComputationClosureSpace, "Computation: Space Closure");

		RDG* rdugSpaceClosed = rdgSpace.calcJoinTransitiveClosure();

		pauseTimer(tcComputationClosureSpace);

		spaceSolutions = new CHadron(calcMultiTupleRelation(*rdugSpaceClosed));

		RDG rdgTime = calculateTime_RDG();
		rdgTime.osHtml = osHtml;

		//RDG* rdugTime = rdgTime.combineConcurrentEdges(ops);
		//rdgTime.getUniqueEdgeBetween(rdgTime.getVertex(0), rdgTime.getVertex(0)).getData().print(true);

		if (tcComputationClosureTime == NULL) { tcComputationClosureTime = new TimeCollector("Computation: Time Closure"); timeCollectors.push_back(tcComputationClosureTime); }
		tcComputationClosureTime->resume();

		RDG* rdugTimeClosed = rdgTime.calcJoinTransitiveClosure();

		tcComputationClosureTime->pause();

		//rdugTimeClosed->getUniqueEdgeBetween(rdugTimeClosed->getVertex(0), rdugTimeClosed->getVertex(0)).getData().print();

		timeSolutions = new CHadron(calcMultiTupleRelation(*rdugTimeClosed));

		//if (osHtml != NULL) *osHtml << "<hr/><hr/><b>Starting Computing Implementation</b><br/>";
		//if (osHtml != NULL) *osHtml << "<h1>Input (Polyhedral representation)</h1><br/>" << toStringHtml() << "<hr/>";

		/* Space-Time
		 * Constant-Affine-...
		 *
		 * Org relation
		 * - For each type of abstraction: {
		 *   - Approximate to space (transitive-reflexive closure)
		 *   - Partition
		 *   - Approximate to time
		 *   - Partition
		 *   }
		 * or
		 * - For each type of abstraction: {
		 *   - Approximate to space
		 *   - Partition
		 *   }
		 * - For each type of abstraction: {
		 *   - Approximate to time
		 *   - Partition
		 *   }
		 */

		// Phase 2: Affine Space Partitioning
/*		for (GraphIterator gI = components->begin(); gI != components->end(); ++gI) { RADGraph& g = **gI;
			H("<b>Component " << compIx++ << "</b><br />");

			newComponents->push_back(g.calcJoinTransitiveClosure());

			levelCDAG.saveToSVGStream(*osHtml);

			H(Html::startTable("Resulting Components", 1));
			vector<RADGraph*> resultingComponents = g.getVertexSetInducedSubGraphs(levelCDAG.getConnectedComponentsData());
			for (GraphIterator hI = resultingComponents.begin(); hI != resultingComponents.end(); ++hI) {
				newComponents->push_back(*hI);
				for (VertexIterator v = (**hI).getVertices().begin(); v != (**hI).getVertices().end(); ++v) {
					if (v != (**hI).getVertices().begin()) H(", "); H(CInteger((*v)->getIndex()).toString());
				}
			} H("</table>");
		}
		delete components; components = newComponents; newComponents = NULL;*/

		// Phase 2: Constant Time Partitioning
		/*newComponents = new vector<RADGraph*>();
		H("<h1>Non-constant Space Partitioning</h1>"); compIx = 0;
		for (GraphIterator gI = components->begin(); gI != components->end(); ++gI) { RADGraph& g = **gI;
			if (osHtml != NULL) *osHtml << "<b>Component " << compIx++ << "</b>";

			CGraph<RADGVertex*, bool> levelACRG = CGraph<RADGVertex*, bool>(); // affine computation relation
			for (VertexIterator v = g.getVertices().begin(); v != g.getVertices().end(); ++v) levelDDG.addVertex(*v);

			H(Html::startTable("Common Data Access Graph", g.getVertices().size()));
			for (VertexIterator v = g.getVertices().begin(); v != g.getVertices().end(); ++v) {  H("<tr>");
				for (VertexIterator w = g.getVertices().begin(); w != g.getVertices().end(); ++w) {  H("<td>");
					vector<RADGEdge*> edges = g.getEdgesBetween(**v, **w);
					if ((edges.size() != 0) && (edges[0] != NULL) && (!edges[0]->getData().isEmpty())) {
						levelDDG.addEdge(levelDDG((*v)->getIndex()), levelDDG((*w)->getIndex()), true);
						H("x");
					} else H("o");
					H("</td>");
				}
				H("</tr>");
			}
			H("</table>");
		}*/
#endif
	}

	/*template <class R>
	CStatement<R>& getYoungestCommonAncestor(const CStatement<R>& sL, const CStatement<R>& sR) {
		vector<CStatement<R>* > ancestorsL;
		vector<CStatement<R>* > ancestorsR;

		CStatement<R>* pL = &sL;
		do {
			ancestorsL.push_back(pL);
			pL = pL->getParent();
		} while (pL != NULL);

		CStatement<R>* pR = &sR;
		do {
			ancestorsL.push_back(pR);
			pR = pR->getParent();
		} while (pR != NULL);

		int cPL = ancestorsL.size() - 1;
		int cPR = ancestorsR.size() - 1;

		while (ancestorsL[cPL] == ancestorsR[cPR]) {
			if ((cPL > 0) && (cPR > 0)) {
				if (ancestorsL[cPL - 1] == ancestorsR[cPR - 1]) {
					cPL--; cPR--;
				}
			}
		}

		return *ancestorsL[cPL];
	}*/

	#define PRE(q) (#q)

	/*template <class R> template <class SetClass>
	CGraph<CStatement<R>*, SetRelation<CFlat<SetClass> > >*
	CComputation<R>::calculateDependenceClosure() {
		typedef CGraph<CStatement<R>*, SetRelation<CFlat<SetClass> > > RADGraph;

		if (osHtml != NULL) *osHtml << "<hr/><hr/><b>Starting Space Partitioning</b><br/>";
		if (osHtml != NULL) *osHtml << "<b>Input (Polyhedral representation)</b><br/>" << toStringHtml() << "<hr/>";

		RADGraph approxRSDG = calculateApproximatedRSDG<SetClass>(false);
		approxRSDG.osHtml = osHtml;

		return approxRSDG.calcJoinTransitiveClosure();
	}*/

	template <class R> void CComputation<R>::analyseCommunicationPartitions() {
	/*	for (unsigned int a = 0; a < variables.size(); a++) { Variable& variable = *variables[a];
			// We will want to reduce the list to unique (statement, transformation) references for performance reasons (later)
			for (int r = 0; r < variable.getReferenceCount(); r++) { const Reference& refR = variable.getReference(r);
				RADGVertex& resultVxV = radg(vxIndices[refR.getStatement().getIndex()]);
				int vxVUDim = refR.getStatement().getFullIterationSpace().getSpaceDimension() - parameterCount;
				const CAffineTransformation<R>& atR = refR.getTransformation();

				for (int s = 0; s <= r; s++) { const Reference& refS = variable.getReference(s);
					RADGVertex& resultVxW = radg(vxIndices[refS.getStatement().getIndex()]);
					int vxWUDim = refS.getStatement().getFullIterationSpace().getSpaceDimension() - parameterCount;
					const CAffineTransformation<R>& atS = refS.getTransformation();

					if (osHtml != NULL) *osHtml << "<tr><td>" << a << "</td><td>" << refR.getStatement().getIndex() << "</td><td>" << refS.getStatement().getIndex() << "</td>";

					// The set relations resulting from the references only
					SetRelation mLR = SetRelation::fromConstraints(vxVUDim, vxWUDim, parameterCount + 1,
						(atR.getSubMatrixColumns(0, affineness) - atS.getSubMatrixColumns(0, affineness))
						<< atR.getSubMatrixColumns(affineness, vxVUDim) << -atS.getSubMatrixColumns(affineness, vxWUDim)
					);
					SetRelation mRL = SetRelation::fromConstraints(vxWUDim, vxVUDim, parameterCount + 1,
						(atS.getSubMatrixColumns(0, affineness) - atR.getSubMatrixColumns(0, affineness))
						<< atS.getSubMatrixColumns(affineness, vxWUDim) << -atR.getSubMatrixColumns(affineness, vxVUDim)
					);

					// Project away the kernel of the iteration space

					radg.addEdge(resultVxV, resultVxW, mLR);
					radg.addEdge(resultVxW, resultVxV, mRL);

					if (osHtml != NULL) *osHtml << "<td>" << mLR.toStringHtml() << "</td><td>" << mRL.toStringHtml() << "</td></tr>";
				}
			}
		}*/
	}

	template <class R> void CComputation<R>::performAffineSetSpacePartitioning() {
		//dependenceFlatClosure = (ModuleRADGraph*) calculateDependenceClosure<CModule<R> >();
	}

	/*template <class R> void CComputation<R>::performAffineLatticeSpacePartitioning() {
		dependenceLatticeClosure = (LatticeRADGraph*) calculateDependenceClosure<CLattice<R> >();

		// Determine distance lattices
		if (osHtml != NULL) *osHtml << "<hr/><b>" << "Calculation of distance lattices:" << "</b>" << "<br/>";
		for (int q = 0; q < dependenceLatticeClosure->getVertexCount(); q++) {
			CRADGLatticeVertex& v = (*dependenceLatticeClosure)(q);
			if (dependenceLatticeClosure->getEdgesBetween(v, v).size() == 1) {
				LatticeRelation l = LatticeRelation(dependenceLatticeClosure->getUniqueEdgeBetween(v, v).getData());
				CLattice<R> kernel = calculateDistanceLattice(l);
				latticeKernels.push_back(kernel);
				if (osHtml != NULL ) *osHtml << "Statement " << q << ":<br/>" << kernel.toStringHtml() << "<br/>";
			}
		}

		// Polyhedral decomposition
		/--decomposedGraph = new CRDMGraph<R>();
		for (int p = 0; p < rdmGraph->getComputationVertexCount(); p++) { CRDMGVertex& vxP = rdmGraph->getComputationVertex(p);
			CRDMGVertex& newVXP = decomposedGraph->addComputationVertex(vxP.getData().polyhedralLatticeDecomposition(*distanceLattices[p]));
			for (int d = 0; d < rdmGraph->getDataVertexCount(); d++) { CRDMGVertex& vxD = rdmGraph->getDataVertex(d);
				if (p == 0) decomposedGraph->addDataVertex(vxD.getData());
				vector<CRDMGEdge* > edges = rdmGraph->getEdgesBetween(vxP, vxD);
				for (unsigned int e = 0; e < edges.size(); e++) { CRDMGEdge& egE = *edges[e];
					CAffineTransformation<R> at = egE.getData().polyhedralLatticeDecomposition(*distanceLattices[p]);
					decomposedGraph->addEdge(newVXP, decomposedGraph->getDataVertex(d), at);
				}
			}
		}--/

		//if (osHtml != NULL) *osHtml << "<hr/><b>Output (Polyhedral representation)</b><br/>" << decomposedGraph->toStringHtml();
	}*/

	template <class R>
	string CComputation<R>::toStringHtml() {
		string result = "<table><tr><td>Statement</td><td>Variable</td><td>Reference</td></tr>";
		for (int s = 0; s < getStatementCount(); s++) { const Statement& statS = getStatement(s);
			for (int r = 0; r < statS.getReferenceCount(); r++) { const Reference& ref = statS.getReference(r);
				result += "<tr><td>" + CInteger(s).toString() + "</td><td>" + CInteger(ref.getVariable().getIndex()).toString() + "</td><td>";
				result += ref.getTransformation().toStringHtml();
				result += "</td></tr>";
			}
		}
		result += "</table>";

		result += "<table><tr><td>Statement</td><td>Iteration Space Constraint Matrix</td></tr>";
		for (int s = 0; s < getStatementCount(); s++) { const Statement& statS = getStatement(s);
			result += "<tr><td>" + CInteger(s).toString() + "</td><td>";
			result += statS.getFullIterationSpace().toStringHtml();
			result += "</td></tr>";
		}
		result += "</table>";

		return result;
	}

	/*
	template <class R> template <class SetClass>
	CGraph<CStatement<R>*, SetRelation<CFlat<SetClass> > >
	CComputation<R>::calculateApproximatedRSDG(bool polyhedral) {
		typedef SetRelation<CFlat<SetClass> > SetRelation;
		typedef CGraph<CStatement<R>*, SetRelation> RADGraph;
		typedef CVertex<CStatement<R>*, SetRelation> RADGVertex;
		typedef typename CComputation<R>::LeafIterator HNLeafIterator;

		if (osHtml != NULL) *osHtml << "Parameters:" << getComputation()->getParameterCount() << "<br/><b>Basic edge relations (each represented as vector space basis):</b><br/><table>"
			<< "<tr><td>Statement 0</td><td>Statement 1</td><td>Variable</td><td>Left statement</td><td>Right statement</td><td>Pairwise relation</td><td>Symmetric relation</td></tr>";

		RADGraph radg = RADGraph();
		radg.osHtml = osHtml;
		for (HNLeafIterator lit = originalExecutionOrder.beginLeaf(); !lit->beyondEnd(); ++lit) {
			radg.addVertex(**lit);
		}

		for (unsigned int a = 0; a < variables.size(); a++) { Variable& variable = *variables[a];
			// We will want to reduce the list to unique (statement, transformation) references for performance reasons (later)
			for (int r = 0; r < variable.getReferenceCount(); r++) { const Reference& refR = variable.getReference(r);
				RADGVertex& resultVxV = radg(refR.getStatement().getIndex());
				for (int s = r; s < variable.getReferenceCount(); s++) { const Reference& refS = variable.getReference(s);
					RADGVertex& resultVxW = radg(refS.getStatement().getIndex());

					if (refR.isWriteReference() || refS.isWriteReference()) {
							if (osHtml != NULL) *osHtml << "<tr><td>" << alit.getIndex() << "</td><td>" << blit.getIndex() << "</td><td>" << a << "</td><td>" << refR.getStatement().getIndex() << "</td><td>" << refS.getStatement().getIndex() << "</td>";

							SetRelation mLR = getSetRelationFromReferencePair<SetClass>(refR, refS, polyhedral);
							SetRelation mRL = getSetRelationFromReferencePair<SetClass>(refS, refR, polyhedral);

							// Project away the kernel of the iteration space

							radg.addEdge(resultVxV, resultVxW, mLR);
							radg.addEdge(resultVxW, resultVxV, mRL);

							if (osHtml != NULL) *osHtml << "<td>" << mLR.toStringHtml() << "</td><td>" << mRL.toStringHtml() << "</td></tr>";
						}
					}
				}
			}
		}

		if (osHtml != NULL) *osHtml << "</table>";

		return radg;
	}*/


	/*template <class R>
	template <class PartitionClass>
	CPartitionRelation<PartitionClass> CComputation<R>::getReferenceSpace(const Reference& fromRef, const Reference& toRef) const {
		CDisjunctiveDomain<Set> fromIt = fromRef.getStatement().getFullIterationSpace();
		CDisjunctiveDomain<Set> toIt = toRef.getStatement().getFullIterationSpace();

		return CPartitionRelation<PartitionClass>::getPairRelationFromCoupledConstraints(
				fromIt.getSpaceDimension() + 1, toIt.getSpaceDimension() + 1,
				fromRef.getTransformation().getConstraintMatrix(),
				toRef.getTransformation().getConstraintMatrix()
		);
	}*/

	template <class R>
	typename CComputation<R>::Set
	CComputation<R>::getConvexRelationDescriptorsFromReferencePair(CReference<R>& refFrom, CReference<R>& refTo) {
		if (tcCRDC == NULL) { tcCRDC = new TimeCollector("ConvexRelation depcone"); timeCollectors.push_back(tcCRDC); }
		tcCRDC->resume();

		CDisjunctiveDomain<Set> dc = getDependenceCone(refFrom, refTo);

		tcCRDC->pause();

		if (tcCRUH == NULL) { tcCRUH = new TimeCollector("ConvexRelation uni hull"); timeCollectors.push_back(tcCRUH); }
		tcCRUH->resume();

		Set dcUH = dc.getUnidirectionalHull();

		tcCRUH->pause();

		return dcUH;
	}


	template <class R>
	typename CComputation<R>::Set
	CComputation<R>::getAffineRelationDescriptorsFromReferencePair(CReference<R>& refFrom, CReference<R>& refTo) {
		if (tcRDC2 == NULL) { tcRDC2 = new TimeCollector("Reflexive Dependence Cone2"); timeCollectors.push_back(tcRDC2); }
		tcRDC2->resume();

		CDisjunctiveDomain<Set> rdc = getReflexiveDependenceCone(refFrom, refTo);

		tcRDC2->pause();

		if (tcRDCBidHull == NULL) { tcRDCBidHull = new TimeCollector("Reflexive Dependence Cone Bid Hull"); timeCollectors.push_back(tcRDCBidHull); }
		tcRDCBidHull->resume();

		Set rdcBH = rdc.getBidirectionalHull();
		//getReflexiveDependenceCone(refFrom, refTo).getBidirectionalHull().print();
		tcRDCBidHull->pause();

		return rdcBH;
	}

	template <class R>
	typename CComputation<R>::DDomain
	CComputation<R>::getRelationDescriptorsFromReferencePair(CReference<R>& fromRef, CReference<R>& toRef) {
		if (tcRDC2 == NULL) { tcRDC2 = new TimeCollector("Reflexive Dependence Cone2"); timeCollectors.push_back(tcRDC2); }
		tcRDC2->resume();

		CDisjunctiveDomain<Set> a = getReferenceSpace(fromRef, toRef);
		CDisjunctiveDomain<Set> b = getPairDomain(fromRef.getStatement(), toRef.getStatement());
		CDisjunctiveDomain<Set> c = getLexiCone(fromRef.getStatement(), toRef.getStatement());
		a.print();
		b.print();
		c.print();
		CDisjunctiveDomain<Set> ab = a && b;
		CDisjunctiveDomain<Set> d = ab && c;

		d.removeEmpties();

		tcRDC2->pause();

		//if (tcRDCBidHull == NULL) { tcRDCBidHull = new TimeCollector("Reflexive Dependence Cone Bid Hull"); timeCollectors.push_back(tcRDCBidHull); }
		//tcRDCBidHull->resume();

		//Set rdcBH = rdc.getBidirectionalHull();
		//getReflexiveDependenceCone(refFrom, refTo).getBidirectionalHull().print();
		//tcRDCBidHull->pause();

		return d;
	}

	template <class R>
	typename CComputation<R>::DDomain
	CComputation<R>::getReflexiveRelationDescriptorsFromReferencePair(CReference<R>& refFrom, CReference<R>& refTo) {
		if (tcRDC2 == NULL) { tcRDC2 = new TimeCollector("Reflexive Dependence Cone2"); timeCollectors.push_back(tcRDC2); }
		tcRDC2->resume();

		CDisjunctiveDomain<Set> rdc = getReflexiveDependenceCone(refFrom, refTo);

		tcRDC2->pause();

		if (tcRDCBidHull == NULL) { tcRDCBidHull = new TimeCollector("Reflexive Dependence Cone Bid Hull"); timeCollectors.push_back(tcRDCBidHull); }
		tcRDCBidHull->resume();

		Set rdcBH = rdc.getBidirectionalHull();
		//getReflexiveDependenceCone(refFrom, refTo).getBidirectionalHull().print();
		tcRDCBidHull->pause();

		return rdcBH;
	}

	template <class R>
	SetRelation<typename CComputation<R>::Set >
	CComputation<R>::getSetRelationFromReferencePair(CReference<R>& refFrom, CReference<R>& refTo, bool polyhedral) {
		return SetRelation<Set>(getParameterCount(), refFrom.getStatement().getIterationSpaceDim(), refTo.getStatement().getIterationSpaceDim(),
				                 polyhedral ? getConvexRelationDescriptorsFromReferencePair(refFrom, refTo)
				                            : getAffineRelationDescriptorsFromReferencePair(refFrom, refTo));
	}

	template <class R>
	SetRelation<typename CComputation<R>::DDomain >
	CComputation<R>::getDisjunctiveSetRelationFromReferencePair(CReference<R>& refFrom, CReference<R>& refTo) {
		return SetRelation<DDomain>(getParameterCount(), refFrom.getStatement().getIterationSpaceDim(), refTo.getStatement().getIterationSpaceDim(),
				                    getRelationDescriptorsFromReferencePair(refFrom, refTo));
	}

	template <class R>
	SetRelation<typename CComputation<R>::DDomain >
	CComputation<R>::getReflexiveDisjunctiveSetRelationFromReferencePair(CReference<R>& refFrom, CReference<R>& refTo) {
		return SetRelation<DDomain>(getParameterCount(), refFrom.getStatement().getIterationSpaceDim(), refTo.getStatement().getIterationSpaceDim(),
				                    getReflexiveRelationDescriptorsFromReferencePair(refFrom, refTo));
	}

	template <class R>
	typename CComputation<R>::Set
	CComputation<R>::getReferenceSpace(Reference& fromRef, Reference& toRef) {
		// XXX: avoid matrix formulation
		CAffineTransformation<R> atF = fromRef.getTransformation();
		CAffineTransformation<R> atT = toRef.getTransformation();
		//atF.print();
		//atT.print();

		int affineness = getParameterCount() + 1;

		//CInteger(fromRef.getStatement().getIterationSpaceDim()).print();

		CMatrix<R> depM
		= (atF.getSubMatrixColumns(0, affineness) - atT.getSubMatrixColumns(0, affineness))
				<< atF.getSubMatrixColumns(affineness, fromRef.getStatement().getIterationSpaceDim())
				<< -atT.getSubMatrixColumns(affineness, toRef.getStatement().getIterationSpaceDim())
		;

		//depM.print();

		typedef typename Set::DeskriptorSetType DST;
		typedef typename Set::DeskriptorType DT;
		DST ds = DST(depM.getColumnCount());
		for (int q = 0; q < depM.getRowCount(); q++) {
			ds.addDeskriptor(DT((Vektor<R>) depM[q], true));
		}

		return Set(C, ds);
	}

	template <class R>
	void CComputation<R>::initCachedDomains() {
		if (!cachedDomainsInitialized) {
			cachedDomainCones = std::vector<std::vector<Set*> >(getStatementCount(), std::vector<Set*>());
			for (int q = 0; q < getStatementCount(); q++) {
				cachedDomainCones[q] = std::vector<Set*>(getStatementCount(), NULL);
			}
		   cachedDomainsInitialized = true;
		}
	}

	template <class R>
	CDisjunctiveDomain<typename CComputation<R>::Set>
	CComputation<R>::getDomainCone(Statement& firstStatement, Statement& secondStatement) {
		initCachedDomains();
		if (cachedDomainCones[firstStatement.getIndex()][secondStatement.getIndex()] != NULL) {
			return *cachedDomainCones[firstStatement.getIndex()][secondStatement.getIndex()];
		}
		//(firstStatement.getFullIterationSpace() ^ secondStatement.getFullIterationSpace()).print();

		//CDisjunctiveDomain<Set> imd = firstStatement.getFullIterationSpace() ^ secondStatement.getFullIterationSpace();
		//imd.print();



		//Deskriptor<R> (*ZV)(int dim, bool bidir) = &(Deskriptor<R>::getZeroDeskriptor);
		//Deskriptor<R> (*UV)(int dim, int unitDim, bool bidir) = &(Deskriptor<R>::getUnitDeskriptor);

		/*DeskriptorSet<Deskriptor<R> > constrs = DeskriptorSet<Deskriptor<R> >(imd.getSpaceDimension() + 1);
		for (int q = 0; q < getParameterCount(); q++) {
			constrs.addDeskriptor(ZV(1, true) << UV(getParameterCount() + firstStatement.getIterationSpaceDim(), q, true)
					                          << -UV(getParameterCount() + secondStatement.getIterationSpaceDim(), q, true));
		}
		//constrs.print();

		imd = imd && Set(Hadron<DeskriptorSet<Deskriptor<R> > >(C, constrs));
		//imd.print();

		vector<int> retainedDimensions = vector<int>();
		for (int q = 0; q < getParameterCount() + 1; q++) retainedDimensions.push_back(q);
		for (int q = 0; q < firstStatement.getIterationSpaceDim(); q++) retainedDimensions.push_back(getParameterCount() + 1 + q);
		for (int q = 0; q < secondStatement.getIterationSpaceDim(); q++) retainedDimensions.push_back(2*getParameterCount() + firstStatement.getIterationSpaceDim() + 1 + q);

		if (tcDomConeProjection == NULL) { tcDomConeProjection = new TimeCollector("Domain Cone Projection"); timeCollectors.push_back(tcDomConeProjection); }
		tcDomConeProjection->resume();

		imd = imd.getProjection(retainedDimensions);*/

		//tcDomConeProjection->pause();
		//imd.print();


		typedef typename Deskriptor<R>::Vector Vektor;

		if (tcFullItSpace == NULL) { tcFullItSpace = new TimeCollector("DC: Full It Space"); timeCollectors.push_back(tcFullItSpace); }
		tcFullItSpace->resume();

		DDomain dsAd = firstStatement.getFullIterationSpace();
		DDomain dsBd = secondStatement.getFullIterationSpace();

		tcFullItSpace->pause();
		if (tcUniHull == NULL) { tcUniHull = new TimeCollector("DC: Unidir hull"); timeCollectors.push_back(tcUniHull); }
		tcUniHull->resume();

		Set dsAh = dsAd.getUnidirectionalHull();
		Set dsBh = dsBd.getUnidirectionalHull();

		tcUniHull->pause();
		if (tcDeskriptors == NULL) { tcDeskriptors = new TimeCollector("DC: desks"); timeCollectors.push_back(tcDeskriptors); }
		tcDeskriptors->resume();

		CDeskriptorSet dsA = dsAh.getDeskriptors(C);
		CDeskriptorSet dsB = dsBh.getDeskriptors(C);

		tcDeskriptors->pause();

		int aff = 1 + getParameterCount();
		int dA = dsA.getDimension() + 1 - aff - 1;
		int dB = dsB.getDimension() + 1 - aff - 1;
		//imd.print();

		Vektor (*ZW)(int dim) = &(Vektor::getZeroVektor);

		CDeskriptorSet result = CDeskriptorSet(aff + dA + dB);
		for (int q = 0; q < dsA.getSize(); q++) {
			Vektor v = dsA[q].getVector().getSubVektor(0, aff) << dsA[q].getVector().getSubVektor(aff, dA) << ZW(dB);
			result.addDeskriptor(Deskriptor<R>(v, dsA[q].isBidirectional()));
		}
		for (int q = 0; q < dsB.getSize(); q++) {
			Vektor v = dsB[q].getVector().getSubVektor(0, aff) << ZW(dA) << dsB[q].getVector().getSubVektor(aff, dB);
			result.addDeskriptor(Deskriptor<R>(v, dsB[q].isBidirectional()));
		}

		CHadron sE = CHadron(C, result);
		Set f = Set(sE);

		cachedDomainCones[firstStatement.getIndex()][secondStatement.getIndex()] = new Set(f);

		return f;

		/*ASSERT_EQ(f, imd.getUnidirectionalHull());

		return imd;*/
	}

	template <class R>
	CDisjunctiveDomain<typename CComputation<R>::Set>
	CComputation<R>::getPairDomain(Statement& firstStatement, Statement& secondStatement) {
		/*initCachedDomains();
		if (cachedDomainCones[firstStatement.getIndex()][secondStatement.getIndex()] != NULL) {
			return *cachedDomainCones[firstStatement.getIndex()][secondStatement.getIndex()];
		}*/

		if (tcFullItSpace == NULL) { tcFullItSpace = new TimeCollector("DC: Full It Space"); timeCollectors.push_back(tcFullItSpace); }
		tcFullItSpace->resume();

		DDomain dsAd = firstStatement.getFullIterationSpace();
		DDomain dsBd = secondStatement.getFullIterationSpace();

		typedef typename Deskriptor<R>::Vector Vektor;
		Vektor (*ZW)(int dim) = &(Vektor::getZeroVektor);

		int aff = 1 + getParameterCount();
		dsAd.print();
		int dA = dsAd.getSpaceDimension() + 1 - aff;
		int dB = dsBd.getSpaceDimension() + 1 - aff;
		DDomain resultDom = DDomain(aff + dA + dB);
		for (int q = 0; q < dsAd.getElementCount(); q++) {
			CDeskriptorSet dsA = dsAd[q].getDeskriptors(C);
			for (int r = 0; r < dsBd.getElementCount(); r++) {
				CDeskriptorSet result = CDeskriptorSet(aff + dA + dB);
				CDeskriptorSet dsB = dsBd[r].getDeskriptors(C);

				for (int q = 0; q < dsA.getSize(); q++) {
					Vektor v = dsA[q].getVector().getSubVektor(0, aff) << dsA[q].getVector().getSubVektor(aff, dA) << ZW(dB);
					result.addDeskriptor(Deskriptor<R>(v, dsA[q].isBidirectional()));
				}
				for (int q = 0; q < dsB.getSize(); q++) {
					Vektor v = dsB[q].getVector().getSubVektor(0, aff) << ZW(dA) << dsB[q].getVector().getSubVektor(aff, dB);
					result.addDeskriptor(Deskriptor<R>(v, dsB[q].isBidirectional()));
				}

				CHadron sE = CHadron(C, result);
				Set f = Set(sE);
				resultDom.addElement(f);
			}
		}

		return resultDom;
	}

	template <class R>
	bool CComputation<R>::doesSyntacticallyPrecede(Statement& firstStatement, Statement& secondStatement) const {
		return transitivelyPrecedes<CComputation<R>, R>(*firstStatement.getLeafExecutionOrderHP(), *secondStatement.getLeafExecutionOrderHP());;
	}

	template <class R>
	CDisjunctiveDomain<typename CComputation<R>::Set>
	CComputation<R>::getLexiCone(Statement& firstStatement, Statement& secondStatement) {
		const ExecutionOrderHP* oncL;
		const ExecutionOrderHP* oncR;

		const ExecutionOrderHP& youngestCommon = getYoungestCommonAndOldestNonCommmonAncestors<CComputation<R>, R>(
				*firstStatement.getLeafExecutionOrderHP(), *secondStatement.getLeafExecutionOrderHP(), &oncL, &oncR);

		//const ExecutionOrderHSP* firstStatementHSP = dynamic_cast<const ExecutionOrderHSP*>(firstStatement.getLeafExecutionOrderHP());
		//const ExecutionOrderHSP* secondStatementHSP = dynamic_cast<const ExecutionOrderHSP*>(secondStatement.getLeafExecutionOrderHP());
		//const ExecutionOrderHSP* youngestCommonHSP = dynamic_cast<const ExecutionOrderHSP*>(&youngestCommon);

		//CInteger(youngestCommon.getIterationSpaceDimension()).print();
		//CInteger(firstStatement.getLeafExecutionOrderHP()->getIterationSpaceDimension()).print();
		//CInteger(secondStatement.getLeafExecutionOrderHP()->getIterationSpaceDimension()).print();
		//if (doesSyntacticallyPrecede(firstStatement, secondStatement)) { printf("yes\n"); } else { printf("no\n"); }
		return getLexiCone(getParameterCount(),
			youngestCommon.getIterationSpaceDimension() - getParameterCount(),
			firstStatement.getIterationSpaceDim(),
			secondStatement.getIterationSpaceDim(),
			doesSyntacticallyPrecede(firstStatement, secondStatement) || (&firstStatement == &secondStatement)
		);
	}

	template <class R>
	CDisjunctiveDomain<typename CComputation<R>::Set>
	CComputation<R>::getLexiCone(int parameterCount, int commonSpaceDim, int firstDim, int secondDim, bool firstSyntPrecedesSecond) {
		typedef DeskriptorSet<Deskriptor<R> > DS;
		CDisjunctiveDomain<Set> result = CDisjunctiveDomain<Set>(firstDim + secondDim + parameterCount + 1);

		DS (*U)(int dim, bool bidir) = &(DS::unitDeskriptorSet);
		DS (*Z)(int count, int dim, bool bidir) = &(DS::zeroDeskriptorSet);
		Deskriptor<R> (*UV)(int dim, int unitDim, bool bidir) = &(Deskriptor<R>::getUnitDeskriptor);
		Deskriptor<R> (*ZV)(int dim, bool bidir) = &(Deskriptor<R>::getZeroDeskriptor);

		// Set = Flat<Hadron<DeskriptorSet<Deskriptor<R> > > >
		for (int q = 0; q < commonSpaceDim; q++) {
			if (q > 0) {
				result.addElement(Set(Hadron<DeskriptorSet<Deskriptor<R> > >(C, DeskriptorSet<Deskriptor<R> >(
					(Z(q, parameterCount + 1, true) <<  U(q, true) << Z(q, firstDim  - q, true) << -U(q, true) << Z(q, secondDim - q, true))
					>> (-UV(parameterCount + 1, 0, false) << -UV(firstDim, q, false) << UV(secondDim, q, false))
					>> UV(parameterCount + 1 + firstDim + secondDim, 0, false)
				))));
			} else {
				if (firstSyntPrecedesSecond) {
					result.addElement(Set(Hadron<DeskriptorSet<Deskriptor<R> > >(C, DeskriptorSet<Deskriptor<R> >(
						DS::fromSingleElement(ZV(parameterCount + 1, false) << -UV(firstDim, q, false) << UV(secondDim, q, false))
					    >> UV(parameterCount + 1 + firstDim + secondDim, 0, false)
					))));
				} else {
					result.addElement(Set(Hadron<DeskriptorSet<Deskriptor<R> > >(C, DeskriptorSet<Deskriptor<R> >(
						DS::fromSingleElement(-UV(parameterCount + 1, 0, false) << -UV(firstDim, q, false) << UV(secondDim, q, false))
					    >> UV(parameterCount + 1 + firstDim + secondDim, 0, false)
					))));
				}
			}
		}

		/*if (firstSyntPrecedesSecond) {
			result.addElement(Set(Hadron<DeskriptorSet<Deskriptor<R> > >(C, DeskriptorSet<Deskriptor<R> >(
				Z(commonSpaceDim, parameterCount +1, true) << U(commonSpaceDim, true) << -U(commonSpaceDim, true)
			))));
		}*/

		return result;
	}

	template <class R>
	CDisjunctiveDomain<typename CComputation<R>::Set> CComputation<R>::getDependenceCone(Reference& fromRef, Reference& toRef) {
		//getReferenceSpace(fromRef, toRef).print();
		//getDomainCone(fromRef.getStatement(), toRef.getStatement()).print();
		//getLexiCone(fromRef.getStatement(), toRef.getStatement()).print();
		//getReflexiveDependenceCone(fromRef, toRef).print();
		//(getReflexiveDependenceCone(fromRef, toRef) && getLexiCone(fromRef.getStatement(), toRef.getStatement())).print();
		/*(getReferenceSpace(fromRef, toRef)
			       && getDomainCone(fromRef.getStatement(), toRef.getStatement())
			       && getLexiCone(fromRef.getStatement(), toRef.getStatement())).print();*/
		//getLexiCone(fromRef.getStatement(), toRef.getStatement()).print();
		return getReflexiveDependenceCone(fromRef, toRef)
		       && getLexiCone(fromRef.getStatement(), toRef.getStatement());
	}


	template <class R>
	CDisjunctiveDomain<typename CComputation<R>::Set> CComputation<R>::getReflexiveDependenceCone(Reference& fromRef, Reference& toRef) {
		if (tcReferenceSpace == NULL) { tcReferenceSpace = new TimeCollector("Reference Space"); timeCollectors.push_back(tcReferenceSpace); }
		tcReferenceSpace->resume();

		CDisjunctiveDomain<Set> refSpace = getReferenceSpace(fromRef, toRef);

		tcReferenceSpace->pause();

		if (tcDomainCone == NULL) { tcDomainCone = new TimeCollector("Domain Cone"); timeCollectors.push_back(tcDomainCone); }
		tcDomainCone->resume();

		CDisjunctiveDomain<Set> domCone = getDomainCone(fromRef.getStatement(), toRef.getStatement());

		tcDomainCone->pause();

		if (tcReflexiveDependenceCone == NULL) { tcReflexiveDependenceCone = new TimeCollector("Reflexive Dependence Cone"); timeCollectors.push_back(tcReflexiveDependenceCone); }
		tcReflexiveDependenceCone->resume();

		CDisjunctiveDomain<Set> result = refSpace && domCone;

		tcReflexiveDependenceCone->pause();
		//printf("RefSpace:\n");
		//getReferenceSpace(fromRef, toRef).print();
		//printf("DomCone:\n");
		//getDomainCone(fromRef.getStatement(), toRef.getStatement()).print();
		/*getLexiCone(fromRef.getStatement(), toRef.getStatement()).print();
		(getReferenceSpace(fromRef, toRef)
			       && getDomainCone(fromRef.getStatement(), toRef.getStatement())
			       && getLexiCone(fromRef.getStatement(), toRef.getStatement())).print();*/
		return result;
	}

	/*template <class R>
	CDisjunctiveDomain<typename CComputation<R>::Set> CComputation<R>::getReflexiveDependenceDomain(Reference& fromRef, Reference& toRef) {
		if (tcReferenceSpace == NULL) { tcReferenceSpace = new TimeCollector("Reference Space"); timeCollectors.push_back(tcReferenceSpace); }
		tcReferenceSpace->resume();

		CDisjunctiveDomain<Set> refSpace = getReferenceSpace(fromRef, toRef);

		tcReferenceSpace->pause();

		if (tcDomainCone == NULL) { tcDomainCone = new TimeCollector("Domain Cone"); timeCollectors.push_back(tcDomainCone); }
		tcDomainCone->resume();

		CDisjunctiveDomain<Set> domCone = getPairDomain(fromRef.getStatement(), toRef.getStatement());

		tcDomainCone->pause();

		if (tcReflexiveDependenceCone == NULL) { tcReflexiveDependenceCone = new TimeCollector("Reflexive Dependence Cone"); timeCollectors.push_back(tcReflexiveDependenceCone); }
		tcReflexiveDependenceCone->resume();

		CDisjunctiveDomain<Set> result = refSpace && domCone;

		tcReflexiveDependenceCone->pause();
		//printf("RefSpace:\n");
		//getReferenceSpace(fromRef, toRef).print();
		//printf("DomCone:\n");
		//getDomainCone(fromRef.getStatement(), toRef.getStatement()).print();
		return result;
	}*/

	template <class R>
	int CComputation<R>::getParameterCount() const {
		if (originalExecutionOrder != NULL) {
			typedef CHierarSetPart<CComputation<R>, R> HierarSetPart;
			HierarSetPart* hierarSetPart = NULL;
			if ((hierarSetPart = dynamic_cast<HierarSetPart* >(originalExecutionOrder)) != NULL) {
				return hierarSetPart->getDomainData().getSpaceDimension();
			} else {
				return 0;
			}
		} else { ASSERT(false); return -1; }  // we probably don't want to be here
	}

}

#endif /*COMPUTATION_H_*/
