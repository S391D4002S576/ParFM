#ifndef CUTE_TEST_H_
#define CUTE_TEST_H_
#include <boost/function.hpp>
#include "cute_listener.h"
// make plain functions as tests more 'cute':
namespace cute {
	class TCuteTest {
	public:
		virtual bool performTest(Listener* L) const = 0;
		
		virtual ~TCuteTest() {};
	};

	class test : public TCuteTest {
	public:
		std::string name()const{ return name_;}
		
		// this shouldn't belong here, but where? needs adaptation for different compilers
		// dependency to demangle is a given, otherwise we have to use macros everywhere
		static std::string demangle(char const *name);
	
		// (real) functor types can (almost) spell their name
		// but a name can also be given explicitely, e.g. for CUTE() macro for simple test functions
		template <typename VoidFunctor>
		test(VoidFunctor const &t, std::string name = demangle(typeid(VoidFunctor).name())) {
			theTest = t;
			name_ = name;
		}
		
		bool performTest(Listener* L) const { 
			const test& t = *this;
			
			try {
				L->start(t);
				t.theTest();
				L->success(t,"OK");
				return true;
			} catch (cute::test_failure const &e) { L->failure(t,e);
			} catch (std::exception const &exc) { L->error(t,test::demangle(exc.what()).c_str());
			} catch (std::string &s) { L->error(t,s.c_str());
			} catch (char const *&cs) { L->error(t,cs);
			} catch (...) { L->error(t,"unknown exception thrown");
			}
			
			return false;
		}
		
		virtual ~test() {};

	private:
		boost::function<void()> theTest;
		std::string name_;
	};
#define CUTE(name) new cute::test((&name),(#name))
	
	// TODO: provide platform independent means of demangling, 
	// or at least support for different compilers
	// this is platform dependant for gnu compilers
#ifdef __GNUG__
#include <cxxabi.h> // __cxa_demangle
	inline std::string test::demangle(char const *name){
		char *toBeFreed = __cxxabiv1::__cxa_demangle(name,0,0,0);
		std::string result(toBeFreed);
		::free(toBeFreed);
		return result;
	}
	
#else
	// this default works reasonably with MSVC71 and 8
	inline std::string test::demangle(char const *name){
		return std::string(name);
	}
#endif
	
}
#endif /*CUTE_TEST_H_*/
