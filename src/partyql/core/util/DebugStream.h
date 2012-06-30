#ifndef DEBUGSTREAM_H
#define DEBUGSTREAM_H

#include <fstream>
#include <string>

	class DebugOutStream {
	public:
		std::string fileName;
		std::string tempData;
		DebugOutStream(std::string iFileName) : fileName(iFileName) {
			tempData = "<html xmlns=\"http://www.w3.org/1999/xhtml\">";
			tempData += "<head>";
			tempData += "<LINK href=\"style.css\" rel=\"stylesheet\" type=\"text/css\">";
			tempData += "<script src=\"fm.js\" type=\"text/javascript\"></script></head>";
			tempData += "<body><font face=\"Bitstream Vera Sans\" size=\"7\">";
		};

		DebugOutStream* operator << (std::string data) {
			tempData += data;
			flush();

			return this;
		}

		void flush() {
			std::ofstream* ofs = new std::ofstream(("examples/" + fileName + ".html").c_str());
			*ofs << tempData;
			for (int q = 0; q < 20; q++) *ofs << "<br>";
			*ofs << "</font></body></html>";
			ofs->close();
		}
	};

extern int svgCount;


#define H(a) { if (osHtml != NULL) *osHtml << (a); }

#endif
