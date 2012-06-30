#ifndef FINITEPARTITION_H_
#define FINITEPARTITION_H_

#include <vector>
using std::vector;

#include <string>
#include <stdio.h>

namespace AlgoTrans {
	class CSet {
	private:
	
	public:
	};
	
	class CFiniteSet : public CSet {
	private:
		int elementCount;
	public:
		CFiniteSet() : elementCount(-1) { };
		CFiniteSet(int iElementCount) : elementCount(iElementCount) { };
		
		int getElementCount() const { return elementCount; }
	};
	
	class CPartition { 
		// Binary ops: coarsen, refine, finerThan, coarserThan, glb, lub, 
	public:
		//virtual const CSet& getSet() const = 0;
	};
	
	class CFinitePartition : public CPartition { // XXX: , public CFiniteSet
		friend bool operator == (const CFinitePartition& a, const CFinitePartition& b);
		friend bool isNoCoarserThan(const CFinitePartition& a, const CFinitePartition& b);
		friend CFinitePartition coarsen(const CFinitePartition& a, const CFinitePartition& b);
		friend CFinitePartition refine(const CFinitePartition& a, const CFinitePartition& b);
	private:
		const CFiniteSet& finiteSet;
		vector<int> kernelFunc;
		
		bool isConsistent() const;
	public:
		CFinitePartition(const CFiniteSet& iFiniteSet, vector<int> iKernelFunc);
		
		const CFiniteSet& getSet() const { return finiteSet; };
		
		int getKernelFuncResult(int element) const;
		
		std::string toString() const;
		void print() const { printf((toString() + "\n").c_str()); };
	};
	
	CFinitePartition coarsen(const CFinitePartition& a, const CFinitePartition& b);
	CFinitePartition operator ||(const CFinitePartition& a, const CFinitePartition& b);
	
	CFinitePartition refine(const CFinitePartition& a, const CFinitePartition& b);
	CFinitePartition operator &&(const CFinitePartition& a, const CFinitePartition& b);
	
	bool isNoCoarserThan(const CFinitePartition& a, const CFinitePartition& b);
	bool operator <= (const CFinitePartition& a, const CFinitePartition& b);
	
	bool operator == (const CFinitePartition& a, const CFinitePartition& b);
	
	class CPartitionRelation {
	public:
		CPartition getPartition(int index);
	};
	
	class CConePartitionRelation {
		
	};
}

#endif
