#include "Html.h"

#include "../scalar/Integer.h"

namespace Html {
    int htmlCurrId = 0;

	std::string wrapper(std::string title, std::string contents, std::string titleBgColor) {
		if (titleBgColor != "") titleBgColor = " bgcolor=\"#" + titleBgColor + "\"";
		return "<table><tr><td colspan=\"2\"" + titleBgColor + ">" + title + "</td></tr>"
		       + "<tr><td>" + contents + "</td></tr></table>";
	}

	std::string wrapInJS(std::string contents, bool floater) {
		htmlCurrId++;
		return "<div id=\"id" + AlgoTrans::CInteger(htmlCurrId).toString() + "\"></div>"
			   + "<script type=\"text/javascript\">document.getElementById(\"id" + AlgoTrans::CInteger(htmlCurrId).toString() + "\").innerHTML = "
				 + (floater ? "Floater(" : "") + "(" + contents + ").toStringHtml()" + (floater ? ")" : "") + ";</script>";
	}

	std::string wrapShort(std::string title, std::string contents, std::string titleBgColor) {
		if (titleBgColor != "") titleBgColor = " bgcolor=\"#" + titleBgColor + "\"";
		return "<table cellspacing=\"0\" cellpadding=\"0\" border=\"1\"><tr><td" + titleBgColor + ">" + title + "</td><td>" + contents + "</td></tr></table>";
	}

	std::string wrapType(std::string title, std::string contents, std::string titleBgColor) {
		return wrapShort(title, contents, titleBgColor);
	}

	std::string startTable(std::string title, int columns) {
		return "<table><tr><td colspan=\"" + AlgoTrans::CInteger(columns).toString() + "\">" + title + "</td></tr>";
	}
}
