#ifndef OSTREAM_LISTENER_H_
#define OSTREAM_LISTENER_H_
#include "cute_listener.h"
#include <iostream>
#include "../src/partyql/core/util/Time.h"

namespace cute {
	// a "root" listener displaying output, may be should be templatized as well
	class ostream_listener : public Listener
	{
		std::ostream &out;
		double tPrev, tPrevSuite;
	public:
		ostream_listener():out(std::cerr){}
		ostream_listener(std::ostream &os):out(os) {} 
		void begin(suite const &t,char const *info){
			out << "beginning: " << info<<std::endl;
			tPrevSuite = rtclock();
		}
		void end(suite const &t, char const *info){
			double tDur = (rtclock() - tPrevSuite);
			out << "ending: " << info << " (" << tDur << "s)" << std::endl;
		}
		void start(test const &t){
			out << "starting: " << t.name()<< std::endl;
			tPrev = rtclock();
		}
		void success(test const &t, char const *msg){
			double tDur = (rtclock() - tPrev);
			out <<  t.name() <<" " << msg << " (" << tDur << "s)" << std::endl;
			printTimeCollectors();
		}
		void failure(test const &t,test_failure const &e){
			out << e.filename << ":" << e.lineno << ": testcase failed: " <<e.reason << " in " << t.name()<< std::endl;
		}
		void error(test const &t, char const *what){
			out << what << " in " << t.name() << std::endl;
		}
	};
}
#endif /*OSTREAM_LISTENER_H_*/
