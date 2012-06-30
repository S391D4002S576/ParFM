#ifndef CUTE_BASE_H_
#define CUTE_BASE_H_
#include <string>
#include <stdio.h>
#include "../src/partyql/core/util/Time.h"

namespace cute{
	struct test_failure {
		std::string reason;
		std::string filename;
		int lineno;
	
		test_failure(std::string const &r,char const *f, int line)
		:reason(r),filename(f),lineno(line)
		{ 	
			printf((r + "\n").c_str());
			printf("Set a breakpoint here.");
			printTimeCollectors();
		}
		
		char const * what() const { return reason.c_str(); }
	};
	
	/*void hypothesisFails(std::string const &r,char const *f, int line)
	{ 	
		printf((r + "\n").c_str());
		printf("Set a breakpoint here.");
	}*/
}
#define ASSERTM(cond,msg) if (!(cond)) throw cute::test_failure((msg),__FILE__,__LINE__)
#define ASSERT(cond) ASSERTM(cond,#cond)

#define ASSERT_EQ(a, b) ASSERTM((a) == (b), "(" + std::string(#a) + "\n== " + (a).toString() + "\n!= " + (b).toString() + "\n== " + std::string(#b) + ")")
#define ASSERT_TYPE(a, b) ASSERTM(dynamic_cast< b >(a) != NULL, "(" + std::string(#a) + " (" + (a).toString() + ") " + "is not of type " + std::string(#b) + ")")

#define FAIL() ASSERTM(false, "FAIL()")
#define FAILM(msg) ASSERTM(false, msg)

#endif /*CUTE_BASE_H_*/
