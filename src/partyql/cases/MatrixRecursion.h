#ifndef MATRIXRECURSION_H_
#define MATRIXRECURSION_H_

namespace AlgoTrans {
	class CMatrixVariable {
	private:
		int index;
		int rowCount, columnCount;
	};
	
	class CMatrixOperator {
	public:
		int getSourceCount() = 0;
		
		
	};
	
	class CMatrixOperation {
	private:
		vector<CMatrixVariable*> sources;
		CMatrixVariable* destination;
		
		CMatrixOperator matrixOperator;
	};

	class CMatrixRecursion {
		
	};
}

#endif /*MATRIXRECURSION_H_*/
