#ifndef CUTE_LISTENER_H_
#define CUTE_LISTENER_H_
#include "cute_base.h"
namespace cute {
	class test;
	class suite;
	
	class Listener {
	public:
		virtual void begin(suite const &s, char const *info) = 0;
		virtual void end(suite const &s, char const *info) = 0;
		virtual void start(test const &t) = 0;
		virtual void success(test const &t,char const *msg) = 0;
		virtual void failure(test const &t,test_failure const &e) = 0;
		virtual void error(test const &t,char const *what) = 0;
		
		virtual ~Listener() {};
	};
}
#endif /* CUTE_LISTENER_H_ */

