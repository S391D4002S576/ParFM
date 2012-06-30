#include "core/Domain.h"
#ifndef HIERARGRAPH_H_
#define HIERARGRAPH_H_

#include "graph/Graph.h"

#include "Declarations.h"

#include <stdarg.h>

#include "StandardOperations.h"
#include "basicmath/scalar/Bool.h"

namespace AlgoTrans {
	template <class UniversalData, class LevelData, class NodeClass> class CHierarNodeLeafIterator;

	template <class UniversalData, class LevelData, class NodeClass>
	class CHierarNode {
		typedef CHierarNode<UniversalData, LevelData, NodeClass> HierarNode;
		typedef typename vector<HierarNode*>::iterator HierarNodeVectorIterator;
		typedef CGraph<HierarNode*, Bool> Graph;
		typedef CVertex<HierarNode*, Bool> Vertex;
		typedef CEdge<HierarNode*, Bool> Edge;

		friend class CHierarNodeLeafIterator<UniversalData, LevelData, NodeClass>;
	protected:
		UniversalData* universalData;
	private:
		LevelData* levelData; // owned by this
		NodeClass* leafData; // none owned by this
		HierarNode* parent;

		Graph subGraph; // each HierarNode is owned by this
		Vertex* vertex; // owned by parent's subGraph

		void setUniversalDataRecursively(UniversalData* newUniversalData);

		void setLeafData(NodeClass* newLeafData) { leafData = newLeafData; if (leafData != NULL) leafData->leafHierarNode = this; }
	public:
		CHierarNode(UniversalData* iUniversalData, LevelData* iLevelData,
				    NodeClass* iLeafData, HierarNode* iParent = NULL)
		: universalData(iUniversalData), levelData(iLevelData), parent(iParent) { setLeafData(iLeafData); };

		CHierarNode(NodeClass* iLeafData)
		: universalData(NULL), levelData(NULL), parent(NULL) { setLeafData(iLeafData); };

		CHierarNode(UniversalData* iUniversalData, HierarNode* iParent = NULL)
		: universalData(iUniversalData), levelData(NULL), leafData(NULL), parent(iParent) { };

		static HierarNode* fromLexicoLinearOrderedSubNodes(UniversalData& iUniversalData, LevelData* iLevelData, vector<HierarNode*> childNodes);
		static HierarNode* fromLexicoLinearOrderedSubNodes(LevelData* iLevelData, vector<HierarNode*> childNodes);
		static HierarNode* fromLexicoLinearOrderedSubNodes(UniversalData& iUniversalData, LevelData* iLevelData, int childCount, ...);
		static HierarNode* fromLexicoLinearOrderedSubNodes(LevelData* iLevelData, int childCount, ...);
		static HierarNode* fromUnorderedSubNodes(vector<HierarNode*> childNodes);

		HierarNode* getParent() const { return parent; }
		int getChildCount() const { return subGraph.getVertexCount(); }
		bool isLeaf() const { return getChildCount() == 0; }

		const NodeClass& getLeafData() const { ASSERT(leafData != NULL); return *leafData; }
		NodeClass& getLeafData() { ASSERT(leafData != NULL); return *leafData; }
		bool hasLevelData() const { return (levelData != NULL); }
		const LevelData& getLevelData() const { ASSERT(hasLevelData()); return *levelData; }
		LevelData& getLevelData() { ASSERT(hasLevelData()); return *levelData; }

		const Vertex& getVertex() const { return *vertex; }

		~CHierarNode() {
			ITT(vector<Vertex* >, v, subGraph.getVertices()) delete (*v)->getData();
			delete levelData;
		}

		typedef CHierarNodeLeafIterator<UniversalData, LevelData, NodeClass> LeafIterator;
		LeafIterator beginLeaf() { return LeafIterator(*this); };
		LeafIterator endLeaf() { return LeafIterator::end(*this); };

		bool precedesAtThisLevel(const HierarNode& a, const HierarNode& b) const;
	};

	template <class UniversalData, class LevelData, class NodeClass>
	const CHierarNode<UniversalData, LevelData, NodeClass>&
	getYoungestCommonAndOldestNonCommmonAncestors(const CHierarNode<UniversalData, LevelData, NodeClass>& sL, const CHierarNode<UniversalData, LevelData, NodeClass>& sR,
			const CHierarNode<UniversalData, LevelData, NodeClass>** oncL, const CHierarNode<UniversalData, LevelData, NodeClass>** oncR);

	template <class UniversalData, class LevelData, class NodeClass>
	class CHierarNodeLeafIterator {
		typedef CHierarNode<UniversalData, LevelData, NodeClass> HierarNode;
		typedef typename vector<CVertex<HierarNode*, Bool>* >::iterator HierarNodeVectorIterator;
	private:
		HierarNode* hierarNode;
		HierarNodeVectorIterator subGraphIterator;
		CHierarNodeLeafIterator* subGraphLeafIterator;
		int index;

		CHierarNodeLeafIterator(HierarNode* iHierarNode,
				HierarNodeVectorIterator iSubGraphIterator,
				CHierarNodeLeafIterator* iSubGraphLeafIterator = NULL, int iIndex = -1)
		: hierarNode(iHierarNode), subGraphIterator(iSubGraphIterator),
		  subGraphLeafIterator(iSubGraphLeafIterator), index(iIndex) { }
	public:
		CHierarNodeLeafIterator() : subGraphLeafIterator(NULL), index(0) {
			subGraphIterator = hierarNode.subGraph.getVertices().begin();
			if (subGraphIterator != hierarNode.subGraph.getVertices().end()) {
				subGraphLeafIterator = new CHierarNodeLeafIterator(**subGraphIterator);
			}
		}

		static CHierarNodeLeafIterator end(HierarNode& hierarNode) {
			return CHierarNodeLeafIterator(hierarNode, hierarNode.subGraphs.end());
		};

		CHierarNodeLeafIterator& operator = (const CHierarNodeLeafIterator& other) {
		    if (this != &other) {
		    	hierarNode = other.hierarNode;
		    	subGraphIterator = other.subGraphIterator;
		    	subGraphLeafIterator = other.subGraphLeafIterator;
		    	index = other.index;
		    }

		    return *this;
		}

		bool operator == (const CHierarNodeLeafIterator& other) const {
			return (subGraphIterator == other.subGraphIterator)
			       && (subGraphLeafIterator == other.subGraphLeafIterator)
			       && (hierarNode == other.hierarNode)
			       && (index == other.index);
		}

		bool operator != (const CHierarNodeLeafIterator& other) const {
			return (subGraphIterator != other.subGraphIterator)
			       || (subGraphLeafIterator != other.subGraphLeafIterator)
			       || (hierarNode != other.hierarNode)
			       || (index != other.index);
		}

		CHierarNodeLeafIterator& operator ++ () {
			if (subGraphIterator != hierarNode.subGraph.getVertices().end()) {
				++subGraphLeafIterator;
				if (subGraphLeafIterator == (*subGraphIterator)->getData()->endLeaf()) {
					++subGraphIterator;
					if (subGraphIterator != hierarNode.end()) {
						subGraphLeafIterator = (*subGraphIterator)->getData()->beginLeaf();
					} else index = -1;
				}
				index++;
			}

			return *this;
		}

		CHierarNodeLeafIterator& operator ++ (int) {
			CHierarNodeLeafIterator tmp = *this;
			++(*this);
			return tmp;
		}

		int getIndex() const { return index; }

		bool beyondEnd() const { return subGraphIterator == hierarNode.end(); }

		NodeClass& operator*() { return (*subGraphLeafIterator)->getData()->getLeafData(); }
		NodeClass* operator->() { return &(*subGraphLeafIterator)->getData()->getLeafData(); }

		CHierarNodeLeafIterator end() { return end(hierarNode); }
	};
}

namespace AlgoTrans {
	template <class UniversalData, class LevelData, class NodeClass>
	const CHierarNode<UniversalData, LevelData, NodeClass>&
	getYoungestCommonAndOldestNonCommmonAncestors(const CHierarNode<UniversalData, LevelData, NodeClass>& sL, const CHierarNode<UniversalData, LevelData, NodeClass>& sR,
			const CHierarNode<UniversalData, LevelData, NodeClass>** oncL, const CHierarNode<UniversalData, LevelData, NodeClass>** oncR) {
		typedef CHierarNode<UniversalData, LevelData, NodeClass> HierarNode;
		vector<const HierarNode* > ancestorsL;
		vector<const HierarNode* > ancestorsR;

		const HierarNode* pL = &sL;
		do { ancestorsL.push_back(pL); pL = pL->getParent(); } while (pL != NULL);

		const HierarNode* pR = &sR;
		do { ancestorsR.push_back(pR); pR = pR->getParent(); } while (pR != NULL);

		int cPL = ancestorsL.size() - 1;
		int cPR = ancestorsR.size() - 1;

		*oncL = ancestorsL[cPL--];
		*oncR = ancestorsR[cPR--];
		while ((*oncL == *oncR) && (cPL >= 0) && (cPR >= 0)) {
			*oncL = ancestorsL[cPL--];
			*oncR = ancestorsR[cPR--];
		}

		return (*oncL == *oncR) ? **oncL : *ancestorsL[cPL + 2];
	}

	template <class UniversalData, class LevelData, class NodeClass>
	bool CHierarNode<UniversalData, LevelData, NodeClass>::precedesAtThisLevel(const HierarNode& a, const HierarNode& b) const {
		typedef PathExistenceOperations<Bool> Operations;

		PathExistenceOperations<Bool> peOps = PathExistenceOperations<Bool>();
		Graph* gFW = subGraph.floydWarshall(peOps);

		std::vector<Edge* > edges = gFW->getEdgesBetween(a.getVertex().getIndex(), b.getVertex().getIndex());

		bool precedes = (edges.size() == 1) && ((bool) edges[0]->getData());

		delete gFW;

		return precedes;
	}

	template <class UniversalData, class LevelData, class NodeClass>
	bool transitivelyPrecedes(const CHierarNode<UniversalData, LevelData, NodeClass>& sL, const CHierarNode<UniversalData, LevelData, NodeClass>& sR) {
		typedef CHierarNode<UniversalData, LevelData, NodeClass> HierarNode;

		const HierarNode* oncL;
		const HierarNode* oncR;

		const HierarNode& ycHN = getYoungestCommonAndOldestNonCommmonAncestors<UniversalData, LevelData, NodeClass>(sL, sR, &oncL, &oncR);

		return ycHN.precedesAtThisLevel(*oncL,*oncR);
	}

	// Ownership of HierarNodes transferred to result
	template <class UniversalData, class LevelData, class NodeClass>
	CHierarNode<UniversalData, LevelData, NodeClass>*
	CHierarNode<UniversalData, LevelData, NodeClass>::fromLexicoLinearOrderedSubNodes(LevelData* iLevelData, vector<HierarNode*> childNodes) {
		HierarNode* result = new HierarNode(NULL, iLevelData, NULL, NULL);

		Vertex* lastV = NULL;
		ITT(vector<HierarNode*>, c, childNodes) {
			(*c)->vertex = &result->subGraph.addVertex(*c);
			if (lastV != NULL) result->subGraph.addEdge(lastV, (*c)->vertex, true);
			lastV = (*c)->vertex;
			(*c)->parent = result;
		}

		return result;
	}

	template <class UniversalData, class LevelData, class NodeClass>
	CHierarNode<UniversalData, LevelData, NodeClass>*
	CHierarNode<UniversalData, LevelData, NodeClass>::fromUnorderedSubNodes(vector<HierarNode*> childNodes) {
		HierarNode* result = new HierarNode(NULL, NULL, NULL, NULL);

		ITT(vector<HierarNode*>, c, childNodes) {
			(*c)->vertex = &result->subGraph.addVertex(*c);
			(*c)->parent = result;
		}

		return result;
	}

	// Ownership of HierarNodes transferred to result
	template <class UniversalData, class LevelData, class NodeClass>
	CHierarNode<UniversalData, LevelData, NodeClass>*
	CHierarNode<UniversalData, LevelData, NodeClass>::fromLexicoLinearOrderedSubNodes(UniversalData& iUniversalData, LevelData* iLevelData, vector<HierarNode*> childNodes) {
		HierarNode* result = fromLexicoLinearOrderedSubNodes(iLevelData, childNodes);

		result->setUniversalDataRecursively(&iUniversalData);

		return result;
	}

	template <class UniversalData, class LevelData, class NodeClass>
	void CHierarNode<UniversalData, LevelData, NodeClass>::setUniversalDataRecursively(UniversalData* newUniversalData) {
		universalData = newUniversalData;
		ITT(vector<Vertex* >, v, subGraph.getVertices()) { (*v)->getData()->setUniversalDataRecursively(newUniversalData); }
	}

	template <class UniversalData, class LevelData, class NodeClass>
	CHierarNode<UniversalData, LevelData, NodeClass>*
	CHierarNode<UniversalData, LevelData, NodeClass>::fromLexicoLinearOrderedSubNodes(UniversalData& iUniversalData, LevelData* iLevelData, int childCount, ...) {
		vector<HierarNode*> childNodes = vector<HierarNode*>();
		va_list argList;
		va_start(argList, childCount);

		while (childCount-- > 0) childNodes.push_back(va_arg(argList, HierarNode*));
		va_end(argList);

		return fromLexicoLinearOrderedSubNodes(iUniversalData, iLevelData, childNodes);
	}

	template <class UniversalData, class LevelData, class NodeClass>
	CHierarNode<UniversalData, LevelData, NodeClass>*
	CHierarNode<UniversalData, LevelData, NodeClass>::fromLexicoLinearOrderedSubNodes(LevelData* iLevelData, int childCount, ...) {
		vector<HierarNode*> childNodes = vector<HierarNode*>();
		va_list argList;
		va_start(argList, childCount);

		while (childCount-- > 0) childNodes.push_back(va_arg(argList, HierarNode*));
		va_end(argList);

		return fromLexicoLinearOrderedSubNodes(iLevelData, childNodes);
	}
}

#endif /*HIERARGRAPH_H_*/
