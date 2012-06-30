#ifndef BOOL_H_
#define BOOL_H_

#include <string>
#include <stdio.h>

namespace AlgoTrans {
	class Bool  {
		friend Bool operator && (const Bool& a, const Bool& b);
		friend Bool operator || (const Bool& a, const Bool& b);
		friend Bool operator + (const Bool& a, const Bool& b); // xor
		friend Bool operator - (const Bool& a, const Bool& b); // xor
		friend bool operator == (const Bool& a, const Bool& b);
	private:
		bool data;

	public:
		Bool() : data(false) { };
		Bool(bool iData) : data(iData) { }
		Bool(int i) { data = (i == 1) ? Bool(true) : Bool(false); }

		operator bool () const { return data; };
		Bool operator ! () const { return Bool(!data); };
		Bool operator - () const { return !(*this); };

		Bool getTransitiveClosure() const { return *this; }

		std::string toString() const { return data ? std::string("1") : std::string("0"); }
		std::string toStringHtml() const { return toString(); }

		void print() const { printf((toString() + "\n").c_str()); };

		Bool& operator += (const Bool& rhs) { data ^= rhs.data; return *this; }

		static std::string getTypeName() { return "Bool"; }
		static std::string getTypeNameHtml() { return "Bool"; }
	};

	Bool operator && (const Bool& a, const Bool& b);
	Bool operator || (const Bool& a, const Bool& b);
	Bool operator + (const Bool& a, const Bool& b);
	Bool operator - (const Bool& a, const Bool& b);
	bool operator == (const Bool& a, const Bool& b);
	bool operator != (const Bool& a, const Bool& b);
}

#endif
