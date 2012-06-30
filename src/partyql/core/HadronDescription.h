#ifndef HADRONDESCRIPTION_H_
#define HADRONDESCRIPTION_H_

namespace AlgoTrans {

//#include "../basicmath/HadronDescription.h"
//typedef HadronDescription Description;

	enum Description { G, C };
	Description operator ! (const Description d);
}

#endif
