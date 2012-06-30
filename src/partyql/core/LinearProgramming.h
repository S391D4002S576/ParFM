#ifndef INTEGERPROGRAMMING_H_
#define INTEGERPROGRAMMING_H_

#include <piplib/piplibMP.h>

#include <string>
using std::string;

#include "scalar/Integer.h"
//#include "Domain.h"
//#include "Descriptor.h"
#include "Vektor.h"
#include "../basicmath/Vektor.h"
#include "scalar/Rational.h"

namespace AlgoTrans {
	//template <class R> class CCone;
	//template <class R> class CModule;
	template <class H> class Flat;
	template <class DS> class Hadron;
	template <class D> class DeskriptorSet;
	template <class R> class Deskriptor;

	template <class R> CVector<typename R::RationalType>* pipVectorToCVector_obs(PipVector * pv);
	template <class R> Vektor<typename R::RationalType>* pipVectorToVektor(PipVector * pv);
	//template <class R> CMatrix<typename R::RationalType>* pipListToCMatrix(PipList * pl);

	template <class R> class PipNewParam {

	};

	template <class R> class CCPipQuast {
		//friend template <class G> CPipQuast<G>* fromPipQuast(PipQuast* iOriginal);

		typedef typename R::RationalType Rational;

		PipNewParam<R>* newParam;

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
		CCPipQuast<R>* nonNegativeCondQuast;
		CCPipQuast<R>* negativeCondQuast;
		CMatrix<Rational>* solution;
		Vektor<Rational>* condition;

		static CCPipQuast<R>* fromPipQuast(PipQuast* iOriginal);

		CCPipQuast() { }

		~CCPipQuast() {
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

	template <class R> class PipProblem {
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

		Flat<Hadron<DeskriptorSet<Deskriptor<R> > > > domConstrFlat;
		Flat<Hadron<DeskriptorSet<Deskriptor<R> > > > domConstrPolyheder;

		Flat<Hadron<DeskriptorSet<Deskriptor<R> > > > paramConstrFlat;
		Flat<Hadron<DeskriptorSet<Deskriptor<R> > > > paramConstrPolyheder;

		PipProblem(int iUnknownCount, int iParameterCount);

		CCPipQuast<R>* solve();
	};

	template <class R> class IpProblem : PipProblem<R> {
	public:
		void setDomainPolyheder(const Flat<Hadron<DeskriptorSet<Deskriptor<R> > > >& newDomainPolyheder) { PipProblem<R>::domConstrPolyheder = newDomainPolyheder; }
		void setDomainEqualities(const Flat<Hadron<DeskriptorSet<Deskriptor<R> > > >& newDomainEqualities) { PipProblem<R>::domConstrFlat = newDomainEqualities; }

		IpProblem(int iUnknownCount) : PipProblem<R>(iUnknownCount, 0) { };

		Vektor<typename R::RationalType>* solve();
	};

	template <class R> class LpProblem : public PipProblem<R> {
	public:
		void setDomainPolyheder(const Flat<Hadron<DeskriptorSet<Deskriptor<R> > > >& newDomainPolyheder) { PipProblem<R>::domConstrPolyheder = newDomainPolyheder; }
		void setDomainEqualities(const Flat<Hadron<DeskriptorSet<Deskriptor<R> > > >& newDomainEqualities) { PipProblem<R>::domConstrFlat = newDomainEqualities; }

		LpProblem(int iUnknownCount) : PipProblem<R>(iUnknownCount, 0) { this->integerSolution = false; };

		Vektor<typename R::RationalType>* solve();
	};

}

#include "IntegerProgramming.cpp"

#endif /*INTEGERPROGRAMMING_H_*/
