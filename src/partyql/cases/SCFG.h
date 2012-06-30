#ifndef SCFG_H_
#define SCFG_H_

#include <vector>

#include "Computation.h"

namespace CFG {
	class Alphabet;
	class Rule;
	class RuleBifurcation;
	class RuleDualEmission;
	class Variable;
	class Grammar;

	class Alphabet {
		int size;

	public:
		Alphabet(int iSize) : size(iSize) { };
	};

	class Rule {
	private:
		bool empty;
	public:
		Rule(bool iEmpty) : empty(iEmpty) { }

		bool isEmpty() const { return empty; }
	};

	class RuleBifurcation : public Rule {
	private:
		Variable* leftVar;
		Variable* rightVar;
	};

	template <class N>
	class RuleDualEmission : public Rule {
	private:
		int leftEmissionLength;
		Variable* centerVar;
		int rightEmissionLength;

		//N getEmissionScore(vector<>) { }
		//std::vector<>
	};

	class Variable {
	private:
		std::vector<Rule*> rules;
		const Grammar& grammar;

		Variable(const Grammar& iGrammar) : grammar(iGrammar) { };
	public:
	};

	class Grammar {
	private:
		std::vector<Variable> variables;
		Alphabet alphabet;
	public:
		Grammar(Alphabet iAlphabet, int numVariables) : alphabet(iAlphabet), variables(std::vector<Variable>(numVariables, Variable(*this))) {  };


	};
}


#endif /* SCFG_H_ */
