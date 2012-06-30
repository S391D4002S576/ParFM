#ifndef POLYHEDRALDOMAIN_CPP_
#define POLYHEDRALDOMAIN_CPP_

#include "PolyhedralDomain.h"
#include "Module.h"
#include "scalar/Interval.h"

#include <list>

namespace AlgoTrans {
#ifdef PPL_LINK
	typedef class Parma_Polyhedra_Library::Powerset<ParmaCPolyhedron> ParmaCPolyhedraPowerSet;
#endif

	template <class R> CPolyhedralDomain<R>::CPolyhedralDomain(const CPolyhedralDomain<R>& iOriginal)
	{
		if (this != &iOriginal) *this = iOriginal;
	}

	template <class R> void CPolyhedralDomain<R>::deletePolyhedra() {
		for (int q = polyhedra.size() - 1; q >= 0; q--) {
			delete polyhedra[q];
		}
		polyhedra.clear();
	}

	template <class R> CPolyhedralDomain<R>& CPolyhedralDomain<R>::operator = (const CPolyhedralDomain<R>& other) {
		if (this != &other) {
			spaceDimension = other.spaceDimension;

			deletePolyhedra();
			for (unsigned int q = 0; q < other.polyhedra.size(); q++) {
				polyhedra.push_back(new CPolyheder<R>(*other.polyhedra[q]));
			}
		}

		return *this;
	}

	template <class R> CInterval<typename R::RationalType> CPolyhedralDomain<R>::asInterval() const {
		typedef typename R::RationalType Rat;
		CInterval<Rat> result = CInterval<Rat>(false, false);

		for (unsigned int q = 0; q < polyhedra.size(); q++) {
			result = cover(result, polyhedra[q]->asInterval());
		}

		return result;
	}

	template <class R> void CPolyhedralDomain<R>::projectFourierMotzkin(int projectionDimension) {
		for (unsigned int q = 0; q < polyhedra.size(); q++) {
			polyhedra[q]->projectFourierMotzkin(projectionDimension);
		}
	}

	template <class R> CPolyhedralDomain<R> CPolyhedralDomain<R>::getProjectionFourierMotzkin(int projectionDimension) const {
		CPolyhedralDomain<R> result = CPolyhedralDomain<R>(getSpaceDimension() - 1);

		for (unsigned int q = 0; q < polyhedra.size(); q++) {
			result.addPolyheder(new CPolyheder<R>(polyhedra[q]->getProjectFourierMotzkin(projectionDimension)));
		}

		return result;
	}

	template <class R> CFlat<CModule<R> > CPolyhedralDomain<R>::getAffineHull() const {
		CFlat<CModule<R> > result = CFlat<CModule<R> >(1, CModule<R>(getSpaceDimension() + 1, G)); // Empty Module


		for (unsigned int q = 0; q < polyhedra.size(); q++) {
			//result.print();
			//polyhedra[q]->getAffineHull().print();
			result = result + polyhedra[q]->getAffineHull();
		}

		return result;
	}

	template <class R> bool CPolyhedralDomain<R>::isEmpty() const {
		for ( int q = 0; q < getPolyhederCount(); q++) {
			if (!getPolyheder(q).isEmpty()) return false;
		}

		return true;
	}

	template <class G> CPolyhedralDomain<G> operator && (const CPolyhedralDomain<G>& a, const CPolyhedralDomain<G>& b) {
		int newSpaceDim = -1;
		if (a.getPolyhederCount() > 0) newSpaceDim = a.getSpaceDimension();
		if (b.getPolyhederCount() > 0) newSpaceDim = b.getSpaceDimension();

		CPolyhedralDomain<G> result = CPolyhedralDomain<G>(newSpaceDim);
		for (int q = 0; q < a.getPolyhederCount(); q++) {
			for (int r = 0; r < b.getPolyhederCount(); r++) {
				result.addPolyheder(a.getPolyheder(q) && b.getPolyheder(r));
			}
		}

		return result;
	}

	template <class G> CPolyhedralDomain<G> operator && (const CPolyheder<G>& a, const CPolyhedralDomain<G>& b) {
		return ((CPolyhedralDomain<G>) a) && b;
	}

	template <class G> CPolyhedralDomain<G> operator && (const CPolyhedralDomain<G>& a, const CPolyheder<G>& b) {
		return a && ((CPolyhedralDomain<G>) b);
	}

	template <class G> CPolyhedralDomain<G> operator || (const CPolyhedralDomain<G>& a, const CPolyhedralDomain<G>& b) {
		CPolyhedralDomain<G> result = CPolyhedralDomain<G>(a.getSpaceDimension());

		for (int q = 0; q < a.getPolyhederCount(); q++) result.addPolyheder(a.getPolyheder(q));
		for (int q = 0; q < b.getPolyhederCount(); q++) result.addPolyheder(b.getPolyheder(q));

		return result;
	}

	template <class R> void CPolyhedralDomain<R>::addPolyheder(const CPolyheder<R>& polyheder) {
		if (getPolyhederCount() == 0) spaceDimension = polyheder.getSpaceDimension();
		ASSERT(polyheder.getSpaceDimension() == spaceDimension);

		polyhedra.push_back(new CPolyheder<R>(polyheder));
	}

	template <class R> CPolyhedralDomain<R> CPolyhedralDomain<R>::getNegative() const {
		CPolyhedralDomain<R> result = CPolyhedralDomain<R>(getSpaceDimension());

		result.addPolyheder(CPolyheder<R>(getSpaceDimension()));

		for (int q = 0; q < getPolyhederCount(); q++) {
			result = result && getPolyheder(q).getNegative();
		}

		return result;
	}

	template <class G> CPolyhedralDomain<G> operator - (const CPolyhedralDomain<G>& a, const CPolyhedralDomain<G>& b) {
		return a && b.getNegative();
	}

	template <class G> bool operator <= (const CPolyhedralDomain<G>& a, const CPolyhedralDomain<G>& b) { // Reflexive SubSetOf
		return (a - b).isEmpty();
	}

	template <class G> bool operator == (const CPolyhedralDomain<G>& a, const CPolyhedralDomain<G>& b) {
		//ASSERT(a.isSpaceCompatibleWith(b));

		return (a <= b) && (b <= a);
	}

	template <class G> bool operator != (const CPolyhedralDomain<G>& a, const CPolyhedralDomain<G>& b) {
		return !(a == b);
	}

	template <class R> void  CPolyhedralDomain<R>::print() const {
		printf((toString() + "\n").c_str());
	}

	template <class R> string CPolyhedralDomain<R>::toString() const {
		string result = "PolyhedralDomain(";

		for (int q = 0; q < getPolyhederCount(); q++) {
			if (q != 0) result += ", ";
			result += getPolyheder(q).toString();
		}

		return result + ")";
	}

	template <class R> string CPolyhedralDomain<R>::toStringHtml() const {
		string result = "<table><tr><td>PolyhedralDomain(</td><td>";

		for (int q = 0; q < getPolyhederCount(); q++) {
			if (q != 0) result += "</td><td>, </td><td>";
			result += getPolyheder(q).toStringHtml();
			result += "</td><td>";
		}

		return result + ")</td></tr></table>";
	}

	template <class R> bool CPolyhedralDomain<R>::isSpaceCompatibleWith(const CPolyhedralDomain<R>& other) const {
		return (getPolyhederCount() == 0) || (other.getPolyhederCount() == 0)
			   || (getSpaceDimension() == other.getSpaceDimension());
	}

	template <class R> CPolyheder<R> CPolyhedralDomain<R>::getConvexHull() const {
		ASSERT(polyhedra.size() > 0);

		//ParmaCPolyhedraPowerSet pResult = *this;
		//pResult.pairwise_reduce();
		bool empty = true;

		CPolyheder<R> result = CPolyheder<R>::source((*polyhedra[0]).getSpaceDimension());
		for (unsigned int q = 0; q < polyhedra.size(); q++) {
			if (!((*polyhedra[q]).isEmpty())) {
				if (empty) {
					result = *polyhedra[q];
					empty = false;
				} else {
					result = result + *polyhedra[q];
				}
			}
		}

		/*CCone<R> dualResult = !CCone<R>::fromConstraints((*polyhedra[0]).getFarkasPolyheder().getConstraintMatrix());
		for (unsigned int q = 0; q < polyhedra.size(); q++) {
			if (!((*polyhedra[q]).isEmpty())) {
				if (empty) {
					dualResult = !CCone<R>::fromConstraints((*polyhedra[q]).getFarkasPolyheder().getConstraintMatrix());
					empty = false;
				} else {
					dualResult = dualResult + !CCone<R>::fromConstraints((*polyhedra[q]).getFarkasPolyheder().getConstraintMatrix());
				}
			}
		}*/

		//result.print();
		//(*pResult.begin()).element().ascii_dump();

		return result;/*
		if (pResult.begin() == pResult.end()) {
			return CPolyheder<R>::source(spaceDimension);
		} else  {
			return CPolyheder<R>::fromParmaCPolyhedron((*pResult.begin()).element());//CPolyheder<R>((!dualResult).getConstraintMatrix());
		}*/
	}

	template <class R> CPolyheder<R> CPolyhedralDomain<R>::getMinimalPositiveHull() const {
		ASSERT(polyhedra.size() > 0);

		bool empty = true;

		CPolyheder<R> result = CPolyheder<R>::source((*polyhedra[0]).getSpaceDimension());
		for (unsigned int q = 0; q < polyhedra.size(); q++) {
			if (!((*polyhedra[q]).isEmpty())) {
				if (empty) {
					result = *polyhedra[q];
					empty = false;
				} else {
					result = result || *polyhedra[q];
				}
			}
		}

		return result;
	}

#ifdef PPL_LINK
	template <class R>
	CPolyhedralDomain<R>::operator ParmaCPolyhedraPowerSet() const {
		ParmaCPolyhedraPowerSet result = ParmaCPolyhedraPowerSet(spaceDimension, Parma_Polyhedra_Library::EMPTY);

		for (unsigned int q = 0; q < polyhedra.size(); q++) {
			//ParmaCPolyhedron pcp = (ParmaCPolyhedron) *polyhedra[q];
			//polyhedra[q]->print();
			//pcp.ascii_dump();

			result.add_disjunct((ParmaCPolyhedron) *polyhedra[q]);
		}

		return result;
	}
#endif
}

#endif
