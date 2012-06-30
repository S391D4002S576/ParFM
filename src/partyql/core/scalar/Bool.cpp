#include "Bool.h"

namespace AlgoTrans {
	Bool operator && (const Bool& a, const Bool& b) { return Bool(((bool) a) && ((bool) b)); };
	Bool operator || (const Bool& a, const Bool& b) { return Bool(((bool) a) || ((bool) b)); };
	Bool operator + (const Bool& a, const Bool& b) { return Bool(((bool) a) ^ ((bool) b)); };
	Bool operator - (const Bool& a, const Bool& b) { return Bool(((bool) a) ^ ((bool) b)); };
	bool operator == (const Bool& a, const Bool& b) { return ((bool) a) == ((bool) b); };
	bool operator != (const Bool& a, const Bool& b) { return !(a == b); };
}
