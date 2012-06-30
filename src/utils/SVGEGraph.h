#ifndef SVGEGRAPH_H_
#define SVGEGRAPH_H_

#include <fstream>
using std::ostream;

#include "SVGElement.h"
#include "../partyql/graph/Graph.h"
using AlgoTrans::CVertex;
using AlgoTrans::CEdge;
using AlgoTrans::CGraph;

#include "SVGECircle.h"
#include "SVGELine.h"
#include "Point.h"

class CSVGEGraph;
class CSVGEGEdgeData;
class CSVGEGVertexData;

class CSVGEGVertexData {
		friend ostream& operator << (ostream &os, CVertex<CSVGEGVertexData, CSVGEGEdgeData>& v);
	public:
		CPoint position;
		CColor fillColor;
		CColor strokeColor;
		int radius;
		int strokeWidth;
		
		CSVGEGVertexData(CPoint iPosition, int iRadius, 
						 CColor iFillColor, CColor iStrokeColor, int iStrokeWidth) 
		: position(iPosition), fillColor(iFillColor), strokeColor(iStrokeColor), radius(iRadius), strokeWidth(iStrokeWidth) { }
		
		virtual ~CSVGEGVertexData() { };
};

class CSVGEGEdgeData {
	friend ostream& operator << (ostream &os, CEdge<CSVGEGVertexData, CSVGEGEdgeData>& e);
	
	public:
		CColor strokeColor;
		int strokeWidth;
		
		CSVGEGEdgeData(CColor iStrokeColor, int iStrokeWidth) : strokeColor(iStrokeColor), strokeWidth(iStrokeWidth) { }
		
		virtual ~CSVGEGEdgeData() { };
};

class CSVGEGraph : public CGraph<CSVGEGVertexData, CSVGEGEdgeData>, public CSVGElement {
	private:
		//int width;
	public:
		//CSVGEGraph(int iWidth) : width(iWidth) { };
		
		void encode(ostream& os) {
			for (int e = 0; e < getEdgeCount(); e++) os << getEdge(e);			
			for (int v = 0; v < getVertexCount(); v++) os << getVertex(v);
		}
		
		virtual ~CSVGEGraph() { };
};

#endif /* SVGEGRAPH_H_ */
