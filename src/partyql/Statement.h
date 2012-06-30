#ifndef STATEMENT_H_
#define STATEMENT_H_

#include "core/Descriptor.h"
#include "core/Domain.h"

#include "graph/Graph.h"

namespace AlgoTrans {
	template <class R> class CComputation;
	template <class R> class CStatement;
	template <class R> class CVariable;
	template <class R> class CReference;

	template <class R> class CModule;
	template <class R> class CAffineTransformation;
	template <class R> class CLattice;

	template <class C, class R> class CHierarPart;
	template <class C, class R> class CHierarSetPart;
	template <class C, class R> class CHierarLeafPart;

	template <class R> class CStatement {
		template <class G> friend void CComputation<G>::registerStatement(CStatement<G>& statement);
		template <class G> friend void CComputation<G>::registerReference(CReference<G>& reference);

		typedef typename CComputation<R>::Set Set;
		typedef CHierarPart<CComputation<R>, R> ExecutionOrderHP;
		friend class CHierarPart<CComputation<R>, R>;
		friend class CHierarLeafPart<CComputation<R>, R>;

		typedef CReference<R> Reference;
	private:
		CDisjunctiveDomain<Set>* cachedFullIterationSpace;

		CComputation<R>& computation;
		vector<Reference*> references;

		int iteratorCount; // For (multi-dimensional) for-loops or similar constructs
		CStatement<R>* parent;
		vector<CStatement<R>*> children;
		CGraph<CStatement<R>*, bool> childDependences;
		CGraph<CStatement<R>*, bool>* childDependencesTransitiveClosure;

		ExecutionOrderHP* leafHierarPart;

	protected:
		int index;

	public:
		CStatement(CComputation<R>& iComputation) : cachedFullIterationSpace(NULL), computation(iComputation), parent(NULL), leafHierarPart(NULL)
		{ computation.registerStatement(*this); }

		// Refs
		int getReferenceCount() const { return references.size(); }
		Reference& getReference(int index) { return *references[index]; }
		void registerReference(Reference& reference);

		// Iterators
		void setIteratorCount(int newIteratorCount) { iteratorCount = newIteratorCount; }
		int getIteratorCount() const { return iteratorCount; }

		void setLeafExecutionOrderHP(ExecutionOrderHP* newLeafExecutionOrderHP) { leafHierarPart = newLeafExecutionOrderHP; };
		ExecutionOrderHP* getLeafExecutionOrderHP() const { return leafHierarPart; };

		// Iteration Space
		int getIterationSpaceDim() { return getFullIterationSpace().getSpaceDimension() - computation.getParameterCount(); };
		CDisjunctiveDomain<Set> getFullIterationSpace();

		int getIndex() const { return index; }
		std::string getName() const { return "s" + CInteger(getIndex()).toString(); }

		~CStatement() { if (cachedFullIterationSpace != NULL) delete cachedFullIterationSpace; }
	};


	/*template <class R>
	class CStatementGroup : public CHierarNode<CComputation<R>, CPolyhedralDomain<R>, CStatement<R> > {
		typedef CStatement<R> Statement;
		typedef CVariable<R> Variable;
		typedef CReference<R> Reference;
	private:
		std::ostream* osHtml;
	public:
		template <class SetClass>
		CGraph<Statement*, CSetRelation<CFlat<SetClass> > >
		calculateApproximatedRSDG(bool polyhedral);

		template <class SetClass> void partitionSpaceConstant();
		template <class SetClass> void partitionSpaceAffine();
		template <class SetClass> void partitionTimeConstant();
		template <class SetClass> void partitionTimeAffine();

		CComputation<R>& getComputation() { return this->universalData; }


	};*/

	/*
	 */

	template <class R> void CStatement<R>::registerReference(CReference<R>& reference) {
		reference.indexInStatement = references.size();

		references.push_back(&reference);
	}

	template <class R>
	CDisjunctiveDomain<typename CStatement<R>::Set> CStatement<R>::getFullIterationSpace() {
		if (tcGetFullIterationSpace == NULL) { tcGetFullIterationSpace = new TimeCollector("Get Full Iteration Space"); timeCollectors.push_back(tcGetFullIterationSpace); }
		tcGetFullIterationSpace->resume();

		if (cachedFullIterationSpace != NULL) return *cachedFullIterationSpace;

		typedef CHierarSetPart<CComputation<R>, R> HierarSetPart;
		//ExecutionOrderHP* currHN  = leafHierarNode;

		//while ((!currHN->hasLevelData()) && (currHN->getParent() != NULL)) currHN = currHN->getParent();


		//int newIts = 0;
		//while (currHN->getParent() != NULL) { currHN = currHN->getParent();
		CDisjunctiveDomain<Set> result = CDisjunctiveDomain<Set>::universe(1);

		bool first = true;
		for (ExecutionOrderHP* currHN = leafHierarPart; currHN != NULL; currHN = currHN->getParent()) {
			HierarSetPart* hierarSetPart = NULL;
			if ((hierarSetPart = dynamic_cast<HierarSetPart* >(currHN)) != NULL) {
				const CDisjunctiveDomain<Set>& currDD = hierarSetPart->getDomainData();
				CDisjunctiveDomain<Set> currDDp = currDD;
				//currDDp.print();
				int newIts = (first) ? 0 : (result.getSpaceDimension() - currDD.getSpaceDimension());
				int dim = (first) ? currDD.getSpaceDimension() : result.getSpaceDimension();
				CDisjunctiveDomain<Set> extPD = CDisjunctiveDomain<Set>(dim);
				ASSERT(first || (extPD.getSpaceDimension() == result.getSpaceDimension()));
				for (int q = 0; q < currDD.getElementCount(); q++) {
					DeskriptorSet<Deskriptor<R> > ds = ((Hadron<DeskriptorSet<Deskriptor<R> > >) currDD[q]).getDeskriptors(C);
					//ds.print();
					DeskriptorSet<Deskriptor<R> > newDS = DeskriptorSet<Deskriptor<R> >(ds.getDimension() + newIts);
					for (int r = 0; r < ds.getSize(); r++) {
						newDS.addDeskriptor(ds[r] << Deskriptor<R>::getZeroDeskriptor(newIts, ds[r].isBidirectional()));
					}
					//Set newSet = Set(Hadron<DeskriptorSet<Deskriptor<R> > >(C, ds << DeskriptorSet<Deskriptor<R> >::zeroDeskriptorSet(ds.getSize(),  newIts, true)));
					//newDS.print();
					//newSet.print();
					extPD.addElement(Set(Hadron<DeskriptorSet<Deskriptor<R> > >(C, newDS)));
					//extPD.print();
				}
				//extPD.print();
				result = (first) ? extPD : (result && extPD);
				//result.print();
				first = false;
			}
		}

		cachedFullIterationSpace = new CDisjunctiveDomain<Set>(result);

		tcGetFullIterationSpace->pause();

		return result;
	}
}

#endif
