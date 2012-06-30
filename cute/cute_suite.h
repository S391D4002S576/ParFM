#ifndef CUTE_SUITE_H_
#define CUTE_SUITE_H_
#include "cute_test.h"
#include "cute_listener.h"
#include <vector>
using std::vector;

namespace cute {
	class suite : public std::vector<TCuteTest*>, public TCuteTest {
	public:
		char const* info;
		
		suite(char const* iInfo) {
			info = iInfo;
		}
		
		bool performTest(Listener* L) const {
			const suite& s = *this;
			L->begin(s, info);
			
			bool result=true;
			for (suite::const_iterator it=s.begin(); it != s.end(); ++it) {
			 	result = (*it)->performTest(L) && result;
			}
			L->end(s, info);
			
			return result;
		}
		
		virtual ~suite() {};
	};
}

#endif /*CUTE_SUITE_H_*/
