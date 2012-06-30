#ifndef SSSVG_H_
#define SSSVG_H_

#include <vector>
using std::vector;

#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <fstream>
using std::ostream;
using std::ofstream;

#include <string>
using std::string;

#include "Point.h"
#include "Color.h"

#include "SVGElement.h"
#include "SVGEElementCollection.h"

#include "SVGELine.h" 
#include "SVGECircle.h"

class CSVG {
	public:
		CSVGElement* rootElement;
		int viewBoxWidth;
		int viewBoxHeight;
		int width;
		int height;
		bool standAlone;
		
		CSVG(int iWidth = 200, int iHeight = 200, bool iStandAlone = false) {
			viewBoxWidth = 50000;
			viewBoxHeight = 50000;
			width = iWidth;
			height = iHeight;
			standAlone = iStandAlone;
		}
		
		CSVG(CSVGElement* iRootElement, int iWidth = 200, int iHeight = 200, bool iStandAlone = false) : rootElement(iRootElement) {
			viewBoxWidth = 50000;
			viewBoxHeight = 50000;
			width = iWidth;
			height = iHeight;
			standAlone = iStandAlone;
		}
		
		void encode(ostream& os) {
			if (standAlone) {
				os << "<?xml version=\"1.0\" standalone=\"no\"?>\n"
				<< "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n";
			}
			
			os
			<< "<svg width=\"" << "100%" //<< "\" heigth=\"" << height 
				<< "\" viewBox=\"0 0 " << viewBoxWidth << " " << viewBoxHeight << "\""
				<< " xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">\n"
			<< rootElement
			<< "</svg>";	
		}
		
		void saveToFile(const string& fileName) {
			ofstream ofs(fileName.c_str());
			
			encode(ofs);
			
			ofs.close();
		}
		
		void saveToPdf(const string& fileName) {
			bool currentStandAlone = standAlone;
			standAlone = true;
			
			string tmpSvg = fileName + ".svg"; 
			saveToFile(tmpSvg);
			
			string tmpEps = fileName + ".eps"; 
			string command = "inkscape --export-eps " + tmpEps + " " + tmpSvg; 
			system(command.c_str());
			remove(tmpSvg.c_str());
						
			string command2 = "epstopdf " + tmpEps + " --outfile=" + fileName; 
			system(command2.c_str());
			remove(tmpEps.c_str());
			
			standAlone = currentStandAlone;
		}
};

#endif /*SSSVG_H_*/
