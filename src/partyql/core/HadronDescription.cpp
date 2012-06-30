#include "HadronDescription.h"

namespace AlgoTrans {

	Description operator ! (const Description d) { return (d == G) ? C : G; }

}
