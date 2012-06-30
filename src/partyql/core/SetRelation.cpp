#ifndef SETRELATION_CPP_
#define SETRELATION_CPP_

#include "SetRelation.h"
#include "scalar/Integer.h"

#include "Hadron.h"

#include "util/Html.h"

namespace AlgoTrans {
	template <class Y> string SetRelation<Y>::toString(bool forceBothDescriptions) {
		return getTypeName() + "(" + CInteger(getPDimension()).toString() + ", " + CInteger(getQDimension()).toString() + ", "
		       + Y::toString(forceBothDescriptions) + ")";
	}

	template <class Y> string SetRelation<Y>::toStringHtml(bool forceBothDescriptions) {
		return Html::wrapInJS(toStringJavaScript(forceBothDescriptions));
/*		if (this == NULL) return "[nil]";
		return Html::wrapper("<table><tr><td>" + getTypeNameHtml()
				             + "</td><td>P-dim: " + CInteger(getPDimension()).toString()
				             + "</td><td>Q-dim: " + CInteger(getQDimension()).toString() + "</td></tr></table>",
		       Y::toStringHtml(forceBothDescriptions), std::string("BBCCAA"));*/
	}

	template <class Y> string SetRelation<Y>::toStringJavaScript(bool forceBothDescriptions) {
		return "new SetRelation(" + CInteger(getPDimension()).toString() + "," + CInteger(getQDimension()).toString() + ","
				+ Y::toStringJavaScript(forceBothDescriptions) + ")";
	}


}

#endif /* SETRELATION_CPP_ */
