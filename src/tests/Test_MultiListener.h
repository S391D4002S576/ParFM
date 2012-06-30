#ifndef TEST_MULTILISTENER_H_
#define TEST_MULTILISTENER_H_

#include "../cute/cute_listener.h"
#include <iostream>
#include <vector>
using std::vector;

namespace cute {
	class Test_MultiListener: public Listener {
	public:
		vector<Listener* > listeners;
		
		Test_MultiListener() {}
			
		void addListener(Listener* L) {
			listeners.push_back(L);
		}
		
		void begin(suite const &t,char const *info){
			for (unsigned int q = 0; q < listeners.size(); q++) listeners[q]->begin(t, info);
		}
		
		void end(suite const &t, char const *info){
			for (unsigned int q = 0; q < listeners.size(); q++) listeners[q]->end(t, info);
		}
		
		void start(test const &t){
			for (unsigned int q = 0; q < listeners.size(); q++) listeners[q]->start(t);
		}
		
		void success(test const &t, char const *msg){
			for (unsigned int q = 0; q < listeners.size(); q++) listeners[q]->success(t, msg);
		}
		
		void failure(test const &t,test_failure const &e){
			for (unsigned int q = 0; q < listeners.size(); q++) listeners[q]->failure(t, e);
		}
		
		void error(test const &t, char const *what){
			for (unsigned int q = 0; q < listeners.size(); q++) listeners[q]->error(t, what);
		}
		
		virtual ~Test_MultiListener() {};
	};
}
#endif /* TEST_MULTILISTENER_H_ */
