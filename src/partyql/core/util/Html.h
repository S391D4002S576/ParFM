#ifndef HTML_H_
#define HTML_H_

#include <string>

namespace Html {
    extern int htmlCurrId;

	std::string wrapInJS(std::string contents, bool floater = true);
	std::string wrapper(std::string title, std::string contents, std::string titleBgColor = "");
	std::string wrapShort(std::string title, std::string contents, std::string titleBgColor = "");
	std::string wrapType(std::string typeName, std::string contents, std::string titleBgColor = "");
	
	std::string startTable(std::string title, int columns);
}

#endif /*HTML_H_*/
