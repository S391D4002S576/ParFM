#ifndef SVGELATTICE_H_
#define SVGELATTICE_H_

#include "../partyql/basicmath/Lattice.h"
using AlgoTrans::CLattice;

#include "../partyql/basicmath/Vektor.h"
using AlgoTrans::CVector;

#include "../partyql/basicmath/scalar/Integer.h"
using AlgoTrans::CInteger;

#include <string>
using std::string;

#include <sstream>
using std::stringstream;

#include "cute.h"
#include "SVGEGraph.h"
#include "SVGELine.h"
#include "SVGEText.h"
#include "ToString.h"

template <class R> class CSVGELattice : public CSVGElement {
	private:
		CLattice<R> lattice;
		CVector<int> lowerBounds;
		CVector<int> upperBounds;
		int width;
		bool flatMode;

		int border;

		int cellDim;
		int xC, yC;

	public:
		string horAxisTitle;
		string verAxisTitle;
		string title;

		bool swapped;

		CSVGELattice(CLattice<R> iLattice, CVector<int> iLowerBounds, CVector<int> iUpperBounds, int iWidth, bool iFlatMode)
		: lattice(iLattice), lowerBounds(iLowerBounds), upperBounds(iUpperBounds), width(iWidth), flatMode(iFlatMode) {
			if (flatMode) {
				ASSERT(lattice.getDimension() - 1 == lowerBounds.getLength());
				ASSERT(lattice.getDimension() - 1 == upperBounds.getLength());

				ASSERT((lattice.getDimension() == 3)
						|| (lattice.getDimension() == 2)); // Only doing these for now
			} else {
				ASSERT(lattice.getDimension() - 1 == lowerBounds.getLength());
				ASSERT(lattice.getDimension() - 1 == upperBounds.getLength());

				ASSERT((lattice.getDimension() == 3)
						|| (lattice.getDimension() == 2)); // Only doing these for now
			}
		};

		CPoint calcGridPos(int x, int y) {
			return CPoint((border + x)*cellDim, (2*yC + border - (border + y))*cellDim);
		}

		CPoint calcCoords(int x, int y) {
			if (lowerBounds.getLength() == 1) {
				return calcGridPos((2*(x - lowerBounds[0]) + 1), 1);
			} else if (lowerBounds.getLength() == 2) {
				return calcGridPos((2*(x - lowerBounds[0]) + 1), 2*(y - lowerBounds[1]) + 1);
			} else return CPoint(0, 0);
		}

		void drawAxis(ostream& os) {
			CColor gridStrokeColor = CColor(0x555555);
			int gridStrokeWidth = (15*(60*cellDim)/100) / 100;
			CColor strokeColor = CColor(0);
			CColor strokeColorMarks = CColor(0);
			int strokeWidth = (25*(60*cellDim)/100) / 100;
			int strokeWidthMarks = (15*(60*cellDim)/100) / 100;
			int markFontSize = (100*cellDim)/100;
			int titleFontSize = (100*cellDim)/100;

			CPoint hatLoc = calcGridPos(2*xC, 0);
			hatLoc.y += 2*cellDim + cellDim/4 + cellDim/4;
			CSVGEText tHorAxisTitle = CSVGEText(hatLoc, strokeColorMarks, horAxisTitle, titleFontSize, "BitStream Vera Sans", "end", "top");

			CPoint vatLoc = calcGridPos(0, 2*yC);
			vatLoc.x -= cellDim + cellDim/2 + cellDim/4 + cellDim/3;
			CSVGEText tVerAxisTitle = CSVGEText(vatLoc, strokeColorMarks, verAxisTitle, titleFontSize, "BitStream Vera Sans", "end", "top");
			//tVerAxisTitle.vertical = true;
			os << &tHorAxisTitle << &tVerAxisTitle;

			CPoint titLoc = calcGridPos(0, 0);
			titLoc.y += cellDim + cellDim/4;
			titLoc.x -= cellDim;
			CSVGEText tTitle = CSVGEText(titLoc, strokeColorMarks, title, titleFontSize, "BitStream Vera Sans", "end", "top");
			os << &tTitle;

			// Marks
			for (int x = lowerBounds[0]; x <= upperBounds[0]; x++) {
				CPoint start = calcGridPos(2*(x - lowerBounds[0]) + 1, 0);
				CPoint end = start;
				end.y += cellDim/4;
				CSVGELine line = CSVGELine(end, start, strokeColorMarks, strokeWidthMarks);

				end.y += cellDim;
				CSVGEText text = CSVGEText(end, strokeColorMarks, intToStr(x), markFontSize, "BitStream Vera Sans", "middle", "top");
				os << &line << &text;

				if (lowerBounds.getLength() == 2) {
					CSVGELine line2 = CSVGELine(start, calcCoords(x, upperBounds[1]),
					                    	gridStrokeColor, gridStrokeWidth);
					os << &line2;
				}
			}

			if (lowerBounds.getLength() == 2) {
				for (int y = lowerBounds[1]; y <= upperBounds[1]; y++) {
					CPoint start = calcGridPos(0, 2*(y - lowerBounds[1]) + 1);
					CPoint end = start;
					end.x -= cellDim/4;
					CSVGELine line = CSVGELine(end, start,
							                    strokeColorMarks, strokeWidthMarks);

					end.x -= cellDim/3;
					CSVGEText text = CSVGEText(end, strokeColorMarks, intToStr(y), markFontSize, "BitStream Vera Sans", "end", "middle");

					CSVGELine line2 = CSVGELine(start, calcCoords(upperBounds[0], y),
							                    gridStrokeColor, gridStrokeWidth);
					os << &line << &text << &line2;
				}
			}

			// Axes
			CPoint zero = calcGridPos(0, 0);
			CSVGELine line = CSVGELine(CPoint(zero.x - strokeWidth/2, zero.y), calcGridPos(2*xC, 0),
					                    strokeColor, strokeWidth);
			os << &line;
			if (lowerBounds.getLength() == 2) {
				CSVGELine line2 = CSVGELine(CPoint(zero.x, zero.y + strokeWidth/2), calcGridPos(0, 2*yC),
					                    strokeColor, strokeWidth);
				os << &line2;
			}
		}

		void encode(ostream& os) {
			border = 3;

			CSVGEGraph g = CSVGEGraph();
			CColor fillColor = CColor(0);
			CColor strokeColor = CColor(0);

			int dim = flatMode ? lattice.getDimension() - 1 : lattice.getDimension();

			xC = upperBounds[0] - lowerBounds[0] + 1;
			cellDim = width / (2*xC + border + 1);
			int radius = ((60 * cellDim) / 100);
			int strokeWidth = ((20 * radius) / 100);

			if (dim == 2) {
				yC = upperBounds[1] - lowerBounds[1] + 1;
				for (int x = lowerBounds[0]; x <= upperBounds[0]; x++) {
					for (int y = lowerBounds[1]; y <= upperBounds[1]; y++) {
						CVector<CInteger> v = CVector<CInteger>();
						if (!swapped) {
							v.appendElement(x);
							v.appendElement(y);
						} else {
							v.appendElement(y);
							v.appendElement(x);
						}
						bool containsPoint = flatMode ? lattice.containsAffinePoint(v) : lattice.containsPoint(v);
						fillColor = containsPoint ? CColor(0xBBAAFF) : CColor(0);
						strokeColor = containsPoint ? CColor(0) : CColor(0);
						int rad = containsPoint ? radius : (radius*2)/7;

						CPoint p = calcCoords(x, y);

						g.addVertex(CSVGEGVertexData(p, rad, fillColor, strokeColor, strokeWidth));
					}
				}
			} else if (dim == 1) {
				yC = 1;
				for (int x = lowerBounds[0]; x <= upperBounds[0]; x++) {
					CVector<CInteger> v = CVector<CInteger>();
					v.appendElement(x);
					bool containsPoint = flatMode ? lattice.containsAffinePoint(v) : lattice.containsPoint(v);
					fillColor = containsPoint ? CColor(0xBBAAFF) : CColor(0);
					strokeColor = containsPoint ? CColor(0) : fillColor;
					int rad = containsPoint ? radius : (radius*2)/7;

					CPoint p = calcCoords(x, 0);

					g.addVertex(CSVGEGVertexData(p, rad, fillColor, strokeColor, strokeWidth));
				}
			}

			drawAxis(os);

			os << &g;
		}

		virtual ~CSVGELattice() { };
};

#endif /* SVGELATTICE_H_ */
