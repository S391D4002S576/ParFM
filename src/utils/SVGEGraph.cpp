#ifndef SVGEGRAPH_CPP_
#define SVGEGRAPH_CPP_

#include "SVGEGraph.h"

ostream& operator << (ostream& os, CEdge<CSVGEGVertexData, CSVGEGEdgeData>& e) {
	CSVGEGVertexData& fromVData = e.getFromVertex().getData();
	CSVGEGVertexData& toVData = e.getToVertex().getData();
	CSVGEGEdgeData& eData = e.getData(); 
	
	CSVGELine(fromVData.position, toVData.position, eData.strokeColor, eData.strokeWidth).encode(os);
	
	return os;
}

ostream& operator << (ostream& os, CVertex<CSVGEGVertexData, CSVGEGEdgeData>& v) {
	const CSVGEGVertexData& data = v.getData();

	CSVGECircle(data.position, data.radius, 
				data.fillColor, data.strokeColor, data.strokeWidth).encode(os);
	
	return os;
}

#endif /* SVGEGRAPH_CPP_ */
