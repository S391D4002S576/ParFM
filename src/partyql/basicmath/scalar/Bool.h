#include "../../core/scalar/Bool.h"

/*
#ifndef BOOL_H_
#define BOOL_H_

#include <string>

namespace AlgoTrans {
	class Bool  {
		friend Bool operator && (const Bool& a, const Bool& b);
		friend Bool operator || (const Bool& a, const Bool& b);
		friend bool operator == (const Bool& a, const Bool& b);
	private:
		bool data;

	public:
		Bool(bool iData) : data(iData) { }

		operator bool () const { return data; };

		std::string toString() const { return data ? std::string("true") : std::string("false"); }
		std::string toStringHtml() const { return toString(); }

		void print() const { printf((toString() + "\n").c_str()); };
	};

	Bool operator && (const Bool& a, const Bool& b) { return Bool(((bool) a) && ((bool) b)); };
	Bool operator || (const Bool& a, const Bool& b) { return Bool(((bool) a) || ((bool) b)); };
	bool operator == (const Bool& a, const Bool& b) { return ((bool) a) == ((bool) b); };
}

#endif
*/
