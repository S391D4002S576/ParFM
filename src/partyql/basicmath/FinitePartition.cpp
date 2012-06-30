#include "FinitePartition.h"

#include "../graph/Graph.h"

#include "../cute/cute.h"
#include "scalar/Integer.h"

namespace AlgoTrans {
	CFinitePartition::CFinitePartition(const CFiniteSet& iFiniteSet, vector<int> iKernelFunc)
	: finiteSet(iFiniteSet), kernelFunc(iKernelFunc) {
		ASSERT(isConsistent());
	}

	int CFinitePartition::getKernelFuncResult(int element) const {
		ASSERT((0 <= element) && (element < finiteSet.getElementCount()));

		return kernelFunc[element];
	}

	bool CFinitePartition::isConsistent() const {
		if (!(finiteSet.getElementCount() == (int) kernelFunc.size())) return false;

		vector<bool> includedCell = vector<bool>(finiteSet.getElementCount(), false);
		int cellCount = 0;
		for (int q = 0; q < finiteSet.getElementCount(); q++) {
			if (!((0 <= kernelFunc[q]) && (kernelFunc[q] < finiteSet.getElementCount()))) return false;

			if (!includedCell[kernelFunc[q]]) {
				includedCell[kernelFunc[q]] = true;
				cellCount++;
			}
		}

		for (int q = 0; q < finiteSet.getElementCount(); q++) {
			if (!((0 <= kernelFunc[q]) && (kernelFunc[q] < cellCount))) return false;
		}

		return true;
	}

	bool operator == (const CFinitePartition& a, const CFinitePartition& b) {
		ASSERT(a.isConsistent() && b.isConsistent());
		ASSERT(a.getSet().getElementCount() == b.getSet().getElementCount());

		vector<int> correspA = vector<int>(a.getSet().getElementCount(), -1);
		vector<int> correspB = vector<int>(a.getSet().getElementCount(), -1);
		for (int q = 0; q < a.getSet().getElementCount(); q++) {
			if (correspA[a.getKernelFuncResult(q)] == -1) {
				if (correspB[b.getKernelFuncResult(q)] != -1) return false;
				correspA[a.getKernelFuncResult(q)] = b.getKernelFuncResult(q);
				correspB[b.getKernelFuncResult(q)] = a.getKernelFuncResult(q);
			} else {
				if (correspA[a.getKernelFuncResult(q)] != b.getKernelFuncResult(q)) return false;
				if (correspB[b.getKernelFuncResult(q)] != a.getKernelFuncResult(q)) return false;
			}
		}

		ASSERT((a <= b) && (b <= a));

		return true;
	}

	CFinitePartition coarsen(const CFinitePartition& a, const CFinitePartition& b) {
		typedef CVertex<int, bool> Vertex;
		ASSERT(a.isConsistent() && b.isConsistent());
		ASSERT(a.getSet().getElementCount() == b.getSet().getElementCount()); // we may actually wan't to check for set equality

		// Setup graph
		CGraph<int, bool> setGraph = CGraph<int, bool>();
		for (int q = 0; q < a.getSet().getElementCount(); q++) setGraph.addVertex(q);

		vector<Vertex* > lastVertices = vector<Vertex* >(a.getSet().getElementCount(), NULL);
		for (int q = 0; q < a.getSet().getElementCount(); q++) {
			if (lastVertices[a.getKernelFuncResult(q)] != NULL) {
				setGraph.addEdge(*lastVertices[a.getKernelFuncResult(q)], setGraph.getVertex(q), true);
			}
			lastVertices[a.getKernelFuncResult(q)] = &setGraph.getVertex(q);
		}

		lastVertices = vector<Vertex* >(b.getSet().getElementCount(), NULL);
		for (int q = 0; q < b.getSet().getElementCount(); q++) {
			if (lastVertices[b.getKernelFuncResult(q)] != NULL) {
				setGraph.addEdge(*lastVertices[b.getKernelFuncResult(q)], setGraph.getVertex(q), true);
			}
			lastVertices[b.getKernelFuncResult(q)] = &setGraph.getVertex(q);
		}

		CConnectedComponentFinder<int, bool> sccf = CConnectedComponentFinder<int, bool>(setGraph);
		vector<vector<Vertex* > > sccs = sccf.findConnectedComponents();

		vector<int> newKernelMap = vector<int>(a.getSet().getElementCount(), -1);
		for (unsigned int scc = 0; scc < sccs.size(); scc++) {
			for (unsigned int v = 0; v < sccs[scc].size(); v++) {
				newKernelMap[sccs[scc][v]->getData()] = scc;
			}
		}

		CFinitePartition result = CFinitePartition(a.getSet(), newKernelMap);

		return result;
	}

	CFinitePartition refine(const CFinitePartition& a, const CFinitePartition& b) {
		typedef CVertex<int, bool> Vertex;
		ASSERT(a.isConsistent() && b.isConsistent());
		ASSERT(a.getSet().getElementCount() == b.getSet().getElementCount()); // we may actually wan't to check for set equality

		// Get a list of included vertices for each cell of A
		vector<int> lastCellVix = vector<int>(a.getSet().getElementCount(), -1);
		vector<int> prevVix = vector<int>(a.getSet().getElementCount(), -1);
		int cellCount = 0;
		for (int q = 0; q < a.getSet().getElementCount(); q++) {
			if (lastCellVix[a.getKernelFuncResult(q)] != -1) {
				prevVix[q] = lastCellVix[a.getKernelFuncResult(q)];
			} else cellCount++;
			lastCellVix[a.getKernelFuncResult(q)] = q;
		}

		// Split each cell into subsets corresponding to cells of b
		vector<int> newKernelMap = vector<int>(a.getSet().getElementCount(), -1);
		int newCellCount = 0;
		for (int c = 0; c < cellCount; c++) {
			vector<int> lastNewVix = vector<int>(a.getSet().getElementCount(), -1);
			for (int v = lastCellVix[c]; v != -1; v = prevVix[v]) {
				if (lastNewVix[b.getKernelFuncResult(v)] == -1) {
					lastNewVix[b.getKernelFuncResult(v)] = newCellCount++;
				}
				newKernelMap[v] = lastNewVix[b.getKernelFuncResult(v)];
			}
		}

		CFinitePartition result =  CFinitePartition(a.getSet(), newKernelMap);

		ASSERT(result.isConsistent());

		return result;
	}

	bool isNoCoarserThan(const CFinitePartition& a, const CFinitePartition& b) {
		typedef CVertex<int, bool> Vertex;
		ASSERT(a.isConsistent() && b.isConsistent());
		ASSERT(a.getSet().getElementCount() == b.getSet().getElementCount()); // we may actually wan't to check for set equality

		// Get a list of included vertices for each cell of A
		vector<int> lastCellVix = vector<int>(a.getSet().getElementCount(), -1);
		vector<int> prevVix = vector<int>(a.getSet().getElementCount(), -1);
		int cellCount = 0;
		for (int q = 0; q < a.getSet().getElementCount(); q++) {
			if (lastCellVix[a.getKernelFuncResult(q)] != -1) {
				prevVix[q] = lastCellVix[a.getKernelFuncResult(q)];
			} else cellCount++;
			lastCellVix[a.getKernelFuncResult(q)] = q;
		}

		// Check that each cell lies entirely in a cell of b
		//int newCellCount = 0;
		for (int c = 0; c < cellCount; c++) {
			if (lastCellVix[c] != -1) {
				int bCell = b.getKernelFuncResult(lastCellVix[c]);
				for (int v = prevVix[lastCellVix[c]]; v != -1; v = prevVix[v]) {
					if (b.getKernelFuncResult(v) != bCell) return false;
				}
			}
		}

		return true;
	}

	std::string CFinitePartition::toString() const {
		std::string result = "[";

		for (int q = 0; q < finiteSet.getElementCount(); q++) {
			if (q != 0) result += ", ";
			result += CInteger(kernelFunc[q]).toString();
		}

		return result + "]";
	}

	CFinitePartition operator ||(const CFinitePartition& a, const CFinitePartition& b) { return coarsen(a, b); };

	CFinitePartition operator &&(const CFinitePartition& a, const CFinitePartition& b) { return refine(a, b); };

	bool operator <=(const CFinitePartition& a, const CFinitePartition& b) { return isNoCoarserThan(a, b); };

}
