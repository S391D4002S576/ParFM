#ifndef HIERARPART_H_
#define HIERARPART_H_

#include <stdarg.h>

#include "Statement.h"

#include "graph/Graph.h"

#include "Declarations.h"
#include "Navigator.h"
#include "Structures.h"

#include "core/scalar/Bool.h"
#include "core/Descriptor.h"

namespace AlgoTrans {
	template <class R> class CReference;
	template <class R> class CComputation;

	template <class UniversalData, class R> class CHierarFinitePart;
	template <class UniversalData, class R> class CHierarSetPart;
	template <class UniversalData, class R> class CHierarScatterPart;

	template <class UniversalData, class R> class CHierarPart {
		friend class CHierarFinitePart<UniversalData, R>;
		friend class CHierarSetPart<UniversalData, R>;
		friend class CHierarScatterPart<UniversalData, R>;

		typedef CHierarPart<UniversalData, R> HierarPart;
	protected:
		typedef CDisjunctiveDomain<Flat<Hadron<DeskriptorSet<Deskriptor<R> > > > > SimpleLevelData;
		UniversalData* universalData;
	private:
		HierarPart* parent;

		virtual void setUniversalDataRecursively(UniversalData* newUniversalData) = 0;
	public:
		CHierarPart(UniversalData* iUniversalData = NULL, HierarPart* iParent = NULL) : universalData(iUniversalData), parent(iParent) { };

		static HierarPart* fromLexicoLinearOrderedSubNodes(UniversalData& iUniversalData, SimpleLevelData iLevelData, vector<HierarPart*> childNodes);
		static HierarPart* fromLexicoLinearOrderedSubNodes(SimpleLevelData iLevelData, vector<HierarPart*> childNodes);
		static HierarPart* fromLexicoLinearOrderedSubNodes(UniversalData& iUniversalData, SimpleLevelData iLevelData, int childCount, ...);
		static HierarPart* fromLexicoLinearOrderedSubNodes(SimpleLevelData iLevelData, int childCount, ...);

		virtual int getChildCount() const = 0;
		virtual bool isLeaf() const { return getChildCount() == 0; };

		HierarPart* getParent() const { return parent; }

		virtual bool precedesAtThisLevel(const HierarPart& a, const HierarPart& b) const { return false; }

		virtual ~CHierarPart() { };

		virtual const CHierarSetPart<UniversalData, R>* getYoungestHierarSetPart() const;
		int getIterationSpaceDimension() const;

		//virtual analyseFromDepth();

		//typedef CHierarPartLeafIterator<UniversalData, R> LeafIterator;
		//LeafIterator beginLeaf() { return LeafIterator(*this); };
		//LeafIterator endLeaf() { return LeafIterator::end(*this); };
	};

	template <class UniversalData, class R> class CHierarFinitePart : public CHierarPart<UniversalData, R> {
		typedef CHierarPart<UniversalData, R> HierarPart;
		typedef CHierarFinitePart<UniversalData, R> HierarFinitePart;
		typedef typename vector<CVertex<HierarPart*, Bool>* >::iterator VertexIt;
	private:
		virtual void setUniversalDataRecursively(UniversalData* newUniversalData);
	public:
		CGraph<HierarPart*, Bool> cellGraph;

		int getChildCount() const { return cellGraph.getVertexCount(); }

		virtual bool precedesAtThisLevel(const HierarPart& a, const HierarPart& b) const;

		static HierarFinitePart* fromLinearOrderedSubNodes(vector<HierarPart*> childNodes);
		static HierarFinitePart* fromUnorderedSubNodes(vector<HierarPart*> childNodes);

		//typedef CHierarLevelFinitePartNavigator<UniversalData, R> HierarLevelNavigator;
		//typedef CHierarFinitePartIterator<UniversalData, R> HierarPartIterator;
		virtual ~CHierarFinitePart() { ITTT(VertexIt, v, cellGraph.getVertices()) delete (*v)->getData(); }
	};

	template <class UniversalData, class R>
	class CHierarSetPart : public CHierarPart<UniversalData, R> {
		typedef CHierarPart<UniversalData, R> HierarPart;
		typedef CHierarSetPart<UniversalData, R> HierarSetPart;
		typedef CDisjunctiveDomain<Flat<Hadron<DeskriptorSet<Deskriptor<R> > > > > SimpleLevelData;
		//typedef CIxMap<Statement*, SimpleLevelData> DomainData;
	private:
		SimpleLevelData domainData;

		HierarPart* childHierarPart;
		bool unordered;

		virtual void setUniversalDataRecursively(UniversalData* newUniversalData) { this->universalData = newUniversalData; childHierarPart->setUniversalDataRecursively(newUniversalData); };
	public:
		CHierarSetPart(SimpleLevelData iDomainData, HierarPart* iChildHierarPart, bool iUnordered)
		: domainData(iDomainData), childHierarPart(iChildHierarPart), unordered(iUnordered) { };

		int getChildCount() const { ASSERT(childHierarPart != NULL); return 1; }

		//virtual const LevelData* getLevelData() const { ASSERT(levelData != NULL); return levelData; }
		//virtual bool hasLevelData() const { return levelData != NULL; }
		const SimpleLevelData& getDomainData() const { return domainData; }
		virtual bool precedesAtThisLevel(const HierarPart& a, const HierarPart& b) const { ASSERT(false); /* only makes sense for finite parts */ }

		static HierarSetPart* fromLexicoOrderedSubNode(SimpleLevelData iDomainData, HierarPart* iChildHierarPart);
		static HierarSetPart* fromUnorderedSubNode(SimpleLevelData iDomainData, HierarPart* iChildHierarPart);

		//typedef CHierarLevelAffinePartNavigator<UniversalData, R> HierarLevelNavigator;
		//typedef CHierarAffinePartIterator<UniversalData, R> HierarPartIterator;
		virtual ~CHierarSetPart() { if (childHierarPart != NULL) delete childHierarPart; };
	};

	template <class UniversalData, class R>
	class CHierarScatterPart : public CHierarPart<UniversalData, R> {
		typedef CHierarPart<UniversalData, R> HierarPart;
		typedef CHierarSetPart<UniversalData, R> HierarSetPart;
		typedef CDisjunctiveDomain<Flat<Hadron<DeskriptorSet<Deskriptor<R> > > > > SimpleLevelData;
		//typedef CIxMap<Statement*, SimpleLevelData> DomainData;
	private:
		typedef Hadron<DeskriptorSet<Deskriptor<R> > > ScatterData;
		ScatterData scatterData;

		HierarPart* childHierarPart;
		bool unordered;

		virtual void setUniversalDataRecursively(UniversalData* newUniversalData) { this->universalData = newUniversalData; childHierarPart->setUniversalDataRecursively(newUniversalData); };
	public:
		CHierarScatterPart(ScatterData iScatterData, HierarPart* iChildHierarPart, bool iUnordered)
		: scatterData(iScatterData), childHierarPart(iChildHierarPart), unordered(iUnordered) { };

		int getChildCount() const { ASSERT(childHierarPart != NULL); return 1; }
		HierarPart* getChildPart() const { return childHierarPart; }
		bool isTimePart() const { return !unordered; }

		//virtual const LevelData* getLevelData() const { ASSERT(levelData != NULL); return levelData; }
		//virtual bool hasLevelData() const { return levelData != NULL; }
		ScatterData& getScatterData() { return scatterData; }
		virtual bool precedesAtThisLevel(const HierarPart& a, const HierarPart& b) const { ASSERT(false); /* only makes sense for finite parts */ }

		static CHierarScatterPart* fromLexicoOrderedSubNode(ScatterData iScatterData, HierarPart* iChildHierarPart);
		static CHierarScatterPart* fromUnorderedSubNode(ScatterData iScatterData, HierarPart* iChildHierarPart);

		//typedef CHierarLevelAffinePartNavigator<UniversalData, R> HierarLevelNavigator;
		//typedef CHierarAffinePartIterator<UniversalData, R> HierarPartIterator;
		virtual ~CHierarScatterPart() { if (childHierarPart != NULL) delete childHierarPart; };
	};

	template <class UniversalData, class R> class CHierarLeafPart : public CHierarPart<UniversalData, R> {
		typedef CHierarPart<UniversalData, R> HierarPart;
		typedef CHierarLeafPart<UniversalData, R> HierarLeafPart;
	private:
		const CStatement<R>& statement;

		virtual void setUniversalDataRecursively(UniversalData* newUniversalData) { this->universalData = newUniversalData; };
	public:
		CHierarLeafPart(CStatement<R>& originalStatement) : statement(originalStatement) { originalStatement.leafHierarPart = this; };

		const CStatement<R>& getStatement() { return statement; }

		int getChildCount() const { return 0; }

		virtual ~CHierarLeafPart() { };
	};

	/*template <class UniversalData, class R>
	class CHierarFinitePartLevelNavigatorControl {
		typedef CHierarPart<UniversalData, R> HierarPart;
		typedef CHierarFinitePart<UniversalData, R> HierarFinitePart;
		//typedef CHierarFinitePartLevelNavigator<UniversalData, R> LevelNavigator;
	private:
		const HierarFinitePart& levelParent;
	public:
		CHierarFinitePartLevelNavigatorControl(const HierarFinitePart& iLevelParent) : levelParent(iLevelParent) { };
	};

	template <class UniversalData, class R>
	class CHierarFinitePartLevelNavigator : public CNavigator<CHierarPart<UniversalData, R> > {
		typedef CHierarPart<UniversalData, R> HierarPart;
		typedef CHierarFinitePart<UniversalData, R> HierarFinitePart;
	private:
	};

	template <class UniversalData, class R>
	class CHierarFinitePartChildIterator {
		typedef CHierarPart<UniversalData, R> HierarPart;
	private:
		typedef typename vector<CVertex<CHierarFinitePart<UniversalData, R>*, Bool>* >::iterator CellIterator;

		const CHierarFinitePart& hierarFinitePart;
		CellIterator cellIt;
		int index;

		CHierarFinitePartChildIterator(const CHierarFinitePart& iHierarFinitePart, CellIterator iCellIt) : hierarFinitePart(iHierarFinitePart), cellIt(iCellIt) {
		}
	public:
		CHierarFinitePartChildIterator(const CHierarFinitePart& iHierarFinitePart) : hierarFinitePart(iHierarFinitePart) {
			cellIt = hierarFinitePart.cellGraph.getVertices().begin();
		}

		static CHierarFinitePartChildIterator begin() {
			return CHierarFinitePartChildIterator(hierarFinitePart);
		}

		static CHierarFinitePartChildIterator end() {
			return CHierarFinitePartChildIterator(hierarFinitePart, hierarFinitePart.cellGraph.getVertices().end());
		}

		CHierarFinitePartChildIterator& operator = (const CHierarFinitePartChildIterator& other) {
		    if (this != &other) {
		    	cellIt = other.cellIt;
		    	hierarFinitePart = other.hierarFinitePart;
		    	index = other.index;
		    }

		    return *this;
		}

		bool operator == (const CHierarFinitePartChildIterator& other) const {
			return (cellIt == other.cellIt) && (hierarFinitePart == other.hierarFinitePart);
		}

		bool operator != (const CHierarFinitePartChildIterator& other) const {
			return (cellIt != other.cellIt) || (hierarFinitePart != other.hierarFinitePart);
		}

		CHierarFinitePartChildIterator& operator ++ () {
			if (cellIt != hierarFinitePart.cellGraph.getVertices().end()) {
				++cellIt;
				index++;
			}

			return *this;
		}

		CHierarFinitePartChildIterator& operator ++ (int) {
			CHierarFinitePartChildIterator tmp = *this;
			++(*this);
			return tmp;
		}

		int getIndex() const { return index; }

		bool beyondEnd() const { return cellIt == hierarFinitePart.cellGraph.getVertices().end(); }

		HierarPart operator*() { return (*cellIt)->getData(); }
		HierarPart* operator->() { return &(*this); }

		CHierarFinitePartChildIterator begin() { return begin(hierarFinitePart); }
		CHierarFinitePartChildIterator end() { return end(hierarFinitePart); }
	};*/


	template <class UniversalData, class R>
	const CHierarSetPart<UniversalData, R>*
	CHierarPart<UniversalData, R>::getYoungestHierarSetPart() const {
		const CHierarPart<UniversalData, R>* currHP = this;

		while (currHP != NULL) {
			const CHierarSetPart<UniversalData, R>* yHSP = dynamic_cast<const CHierarSetPart<UniversalData, R>* >(currHP);
			if (yHSP != NULL) { return yHSP; }
			else { currHP = currHP->parent; }
		}

		return NULL;
	}

	template <class UniversalData, class R>
	int CHierarPart<UniversalData, R>::getIterationSpaceDimension() const {
		const CHierarSetPart<UniversalData, R>* yHSP = getYoungestHierarSetPart();

		if (yHSP != NULL) { return yHSP->getDomainData().getSpaceDimension(); }
		else { return -1; }
	}

	// Ownership of HierarParts transferred to result
	template <class UniversalData, class R>
	CHierarPart<UniversalData, R>*
	CHierarPart<UniversalData, R>::fromLexicoLinearOrderedSubNodes(UniversalData& iUniversalData, SimpleLevelData iLevelData, vector<HierarPart*> childNodes) {
		CHierarPart<UniversalData, R>* result = fromLexicoLinearOrderedSubNodes(iLevelData, childNodes);

		result->setUniversalDataRecursively(&iUniversalData);

		return result;
	}

	// Ownership of HierarParts transferred to result
	template <class UniversalData, class R>
	CHierarPart<UniversalData, R>*
	CHierarPart<UniversalData, R>::fromLexicoLinearOrderedSubNodes(SimpleLevelData iLevelData, vector<HierarPart*> childNodes) {
		return CHierarSetPart<UniversalData, R>::fromLexicoOrderedSubNode(iLevelData,
				CHierarFinitePart<UniversalData, R>::fromLinearOrderedSubNodes(childNodes));
	}

	// Ownership of HierarParts transferred to result
	template <class UniversalData, class R>
	CHierarFinitePart<UniversalData, R>*
	CHierarFinitePart<UniversalData, R>::fromLinearOrderedSubNodes(vector<HierarPart*> childNodes) {
		HierarFinitePart* result = new HierarFinitePart();

		CVertex<HierarPart*, Bool>* lastV = NULL;
		ITT(vector<HierarPart*>, c, childNodes) {
			CVertex<HierarPart*, Bool>* newV =
			/*(*c)->vertex =*/ &result->cellGraph.addVertex(*c);
			if (lastV != NULL) result->cellGraph.addEdge(lastV, newV, true);
			lastV = newV;
			(*c)->parent = result;
		}

		return result;
	}

	template <class UniversalData, class R>
	void CHierarFinitePart<UniversalData, R>::setUniversalDataRecursively(UniversalData* newUniversalData) {
		this->universalData = newUniversalData;
		for (int q = 0; q < cellGraph.getVertexCount(); q++) {
			cellGraph(q).getData()->setUniversalDataRecursively(newUniversalData);
		}
	}

	// Ownership of HierarParts transferred to result
	template <class UniversalData, class R>
	CHierarFinitePart<UniversalData, R>*
	CHierarFinitePart<UniversalData, R>::fromUnorderedSubNodes(vector<HierarPart*> childNodes) {
		HierarFinitePart* result = new HierarFinitePart();

		ITT(vector<HierarPart*>, c, childNodes) {
			/*(*c)->vertex =*/ &result->cellGraph.addVertex(*c); // XXX: vertex required?
			(*c)->parent = result;
		}

		return result;
	}

	template <class UniversalData, class R>
	CHierarSetPart<UniversalData, R>*
	CHierarSetPart<UniversalData, R>::fromUnorderedSubNode(SimpleLevelData iDomainData, HierarPart* iChildHierarPart) {
		CHierarSetPart<UniversalData, R>* result = new CHierarSetPart<UniversalData, R>(iDomainData, iChildHierarPart, true);

		iChildHierarPart->parent = result;

		return result;
	}

	template <class UniversalData, class R>
	CHierarSetPart<UniversalData, R>*
	CHierarSetPart<UniversalData, R>::fromLexicoOrderedSubNode(SimpleLevelData iDomainData, HierarPart* iChildHierarPart) {
		CHierarSetPart<UniversalData, R>* result = new CHierarSetPart<UniversalData, R>(iDomainData, iChildHierarPart, false);

		iChildHierarPart->parent = result;

		return result;
	}

	template <class UniversalData, class R>
	CHierarScatterPart<UniversalData, R>*
	CHierarScatterPart<UniversalData, R>::fromUnorderedSubNode(ScatterData iScatterData, HierarPart* iChildHierarPart) {
		CHierarScatterPart<UniversalData, R>* result = new CHierarScatterPart<UniversalData, R>(iScatterData, iChildHierarPart, true);

		iChildHierarPart->parent = result;

		return result;
	}

	template <class UniversalData, class R>
	CHierarScatterPart<UniversalData, R>*
	CHierarScatterPart<UniversalData, R>::fromLexicoOrderedSubNode(ScatterData iScatterData, HierarPart* iChildHierarPart) {
		CHierarScatterPart<UniversalData, R>* result = new CHierarScatterPart<UniversalData, R>(iScatterData, iChildHierarPart, false);

		iChildHierarPart->parent = result;

		return result;
	}

	template <class UniversalData, class R>
	CHierarPart<UniversalData, R>*
	CHierarPart<UniversalData, R>::fromLexicoLinearOrderedSubNodes(UniversalData& iUniversalData, SimpleLevelData iLevelData, int childCount, ...) {
		vector<HierarPart*> childNodes = vector<HierarPart*>();
		va_list argList;
		va_start(argList, childCount);

		while (childCount-- > 0) childNodes.push_back(va_arg(argList, HierarPart*));
		va_end(argList);

		return fromLexicoLinearOrderedSubNodes(iUniversalData, iLevelData, childNodes);
	}

	template <class UniversalData, class R>
	CHierarPart<UniversalData, R>*
	CHierarPart<UniversalData, R>::fromLexicoLinearOrderedSubNodes(SimpleLevelData iLevelData, int childCount, ...) {
		vector<HierarPart*> childNodes = vector<HierarPart*>();
		va_list argList;
		va_start(argList, childCount);

		while (childCount-- > 0) childNodes.push_back(va_arg(argList, HierarPart*));
		va_end(argList);

		return fromLexicoLinearOrderedSubNodes(iLevelData, childNodes);
	}

	template <class UniversalData, class R>
	const CHierarPart<UniversalData, R>&
	getYoungestCommonAndOldestNonCommmonAncestors(const CHierarPart<UniversalData, R>& sL, const CHierarPart<UniversalData, R>& sR,
			const CHierarPart<UniversalData, R>** oncL, const CHierarPart<UniversalData, R>** oncR) {
		typedef CHierarPart<UniversalData, R> HierarPart;
		vector<const HierarPart* > ancestorsL;
		vector<const HierarPart* > ancestorsR;

		const HierarPart* pL = &sL;
		do { ancestorsL.push_back(pL); pL = pL->getParent(); } while (pL != NULL);

		const HierarPart* pR = &sR;
		do { ancestorsR.push_back(pR); pR = pR->getParent(); } while (pR != NULL);

		int cPL = ancestorsL.size() - 1;
		int cPR = ancestorsR.size() - 1;

		*oncL = ancestorsL[cPL--];
		*oncR = ancestorsR[cPR--];
		while ((*oncL == *oncR) && (cPL >= 0) && (cPR >= 0)) {
			*oncL = ancestorsL[cPL--];
			*oncR = ancestorsR[cPR--];
		}

		if (*oncL == *oncR) return **oncL;

		return *ancestorsL[cPL + 2];
	}

	class PathExistenceOperations {
	public:
		static Bool calcTransition(Bool& a, Bool& b) { return a && b; }
		static Bool calcCombination(Bool& a, Bool& b) { return a || b; }
	};


	template <class UniversalData, class R>
	bool CHierarFinitePart<UniversalData, R>::precedesAtThisLevel(const HierarPart& a, const HierarPart& b) const {
		PathExistenceOperations peOps = PathExistenceOperations();
		CGraph<HierarPart*, Bool>* gFW = cellGraph.floydWarshall(peOps);

		// First find vertex indices
		int aI = -1, bI = -1;
		for (unsigned int q = 0; q < gFW->getVertices().size(); q++) { const CVertex<HierarPart*, Bool>& v = gFW->getVertex(q);
			if (v.getData() == &a) { ASSERT(aI == -1); aI = v.getIndex(); }
			if (v.getData() == &b) { ASSERT(bI == -1); bI = v.getIndex(); }
		} ASSERT((aI != -1) && (bI != -1));

		std::vector<CEdge<HierarPart*, Bool>* > edges = gFW->getEdgesBetween(aI, bI);

		ASSERT(edges.size() <= 1);
		bool precedes = (edges.size() == 1) && ((bool) edges[0]->getData());

		delete gFW;

		return precedes;
	}

	template <class UniversalData, class R>
	bool transitivelyPrecedes(const CHierarPart<UniversalData, R>& sL, const CHierarPart<UniversalData, R>& sR) {
		typedef CHierarPart<UniversalData, R> HierarPart;

		const HierarPart* oncL;
		const HierarPart* oncR;

		const HierarPart& ycHN = getYoungestCommonAndOldestNonCommmonAncestors<UniversalData, R>(sL, sR, &oncL, &oncR);

		return ycHN.precedesAtThisLevel(*oncL,*oncR) || (&sL == &sR);
	}

}

#endif /*HIERARGRAPH_H_*/
