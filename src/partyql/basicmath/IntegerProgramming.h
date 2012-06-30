#ifndef INTEGERPROGRAMMING2_H_
#define INTEGERPROGRAMMING2_H_

#include <piplib/piplibMP.h>
//#include <piplib/type.h>

#include <string>
using std::string;

#include "scalar/Integer.h"
#include "Matrix.h"
#include "Flat.h"
#include "scalar/Rational.h"

namespace AlgoTrans {
	template <class R> class CCone;
	template <class R> class CModule;

	template <class R> CVector<typename R::RationalType>* pipVectorToCVector(PipVector * pv);
	template <class R> CMatrix<typename R::RationalType>* pipListToCMatrix(PipList * pl);

	template <class R> class CPipNewParam {

	};

	template <class R> class CPipQuast {
		//friend template <class G> CPipQuast<G>* fromPipQuast(PipQuast* iOriginal);

		typedef typename R::RationalType Rational;

		CPipNewParam<R>* newParam;

	public:
		CMatrix<Rational>* getBigParmSolution(int bpIndex) {
			if (condition != NULL) {
				if ((*condition)[bpIndex] > Rational::getZero()) {
					return nonNegativeCondQuast->getBigParmSolution(bpIndex);
				} else {
					return negativeCondQuast->getBigParmSolution(bpIndex);
				}
			} else return solution;
		}
	public:
		CPipQuast<R>* nonNegativeCondQuast;
		CPipQuast<R>* negativeCondQuast;
		CMatrix<Rational>* solution;
		CVector<Rational>* condition;

		static CPipQuast<R>* fromPipQuast(PipQuast* iOriginal);

		CPipQuast() { }

		~CPipQuast() {
			delete nonNegativeCondQuast;
			delete negativeCondQuast;
			delete solution;
			delete condition;
		}

		string toString(string indent = "") const {
			string result = indent + "(";
	    	/*if (newParams.length > 0) {
	    		sb.append("\n");
	    		for (int q = 0; q < newParams.length; q++) {
	    			sb.append(indentor);
	    			sb.append(newParams[q].convertToString());
	    			sb.append("\n");
	    		}
	    	}*/

		    if (condition != NULL) {
	    		result += "if (" + condition->toString() + ")\n";
	    		result += nonNegativeCondQuast->toString(indent + "  ");
	    		result += negativeCondQuast->toString(indent + "  ");
	    		result += indent;
	    	} else {
	    		if (solution == NULL) { result += "[Solution == NULL ???]"; }
	    		else result += solution->toString();
	    	}

		    return result + ")\n";
		}

		void print() { printf((toString()).c_str()); }
	};

	template <class R> class CPipProblem {
	private:
		PipMatrix* getPipConstrMatrixFromAffineSetAndPolyheder(int colCount, const CMatrix<R>& eqs, const CMatrix<R>& ineqs);

		int unknownCount, parameterCount;
	public:
		// Options
		int bigParameter;
		bool integerSolution; // false = rational
		bool simplify;
		bool deepestCut;
		bool positiveParameters;
		bool positiveUnknowns;
		bool maximize;

		CFlat<CModule<R> > domConstrFlat;
		CFlat<CCone<R> > domConstrPolyheder;

		CFlat<CModule<R> > paramConstrFlat;
		CFlat<CCone<R> > paramConstrPolyheder;

		CPipProblem(int iUnknownCount, int iParameterCount);

		CPipQuast<R>* solve();
	};

	template <class R> class CIpProblem : CPipProblem<R> {
	public:
		void setDomainPolyheder(const CFlat<CCone<R> >& newDomainPolyheder) { CPipProblem<R>::domConstrPolyheder = newDomainPolyheder; }

		CIpProblem(int iUnknownCount) : CPipProblem<R>(iUnknownCount, 0) { };

		CVector<typename R::RationalType>* solve();
	};

	template <class R> class CLpProblem : public CPipProblem<R> {
	public:
		void setDomainPolyheder(const CFlat<CCone<R> >& newDomainPolyheder) { CPipProblem<R>::domConstrPolyheder = newDomainPolyheder; }

		CLpProblem(int iUnknownCount) : CPipProblem<R>(iUnknownCount, 0) { this->integerSolution = false; };

		CVector<typename R::RationalType>* solve();
	};

}

#include "IntegerProgramming.cpp"

#endif /*INTEGERPROGRAMMING_H_*/
