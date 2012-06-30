#ifndef CUTE_COUNTING_LISTENER_H_
#define CUTE_COUNTING_LISTENER_H_
#include "cute_listener.h"
namespace cute{
	class counting_listener : public Listener {
	public:
		counting_listener()
		:Listener()
		,numberOfTests(0),successfulTests(0),failedTests(0),errors(0),numberOfSuites(0){}
	
		counting_listener(Listener const &s)
		:Listener(s)
		,numberOfTests(0),successfulTests(0),failedTests(0),errors(0),numberOfSuites(0){}

		void begin(suite const &s, char const *info){
			++numberOfSuites;
		}
		void start(test const &t){
			++numberOfTests;
		}
		void success(test const &t,char const *msg){
			++successfulTests;
		}
		void failure(test const &t,test_failure const &e){
			++failedTests;
		}
		void error(test const &t,char const *what){
			++errors;
		}
		void end(suite const &s, char const *info){
		}
		
		int numberOfTests;
		int successfulTests;
		int failedTests;
		int errors;
		int numberOfSuites;
		
		virtual ~counting_listener() {};
	};
}
#endif /*CUTE_COUNTING_LISTENER_H_*/
