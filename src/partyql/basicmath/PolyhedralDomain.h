#ifndef POLYHEDRALDOMAIN_H_
#define POLYHEDRALDOMAIN_H_

#include <vector>
using std::vector;

#include <list>
using std::list;

#include <string>
using std::string;

//#include <ppl.hh>
#include "Polyheder.h"

namespace AlgoTrans {
	template <class R> class CModule;
	template <class R> class CInterval;
	//template <class R> class CPolyheder;
	template <class R> class CLattice;
	template <class R> class CFlat;
	template <class R> class CCone;
	template <class R> class CPolyheder;

	template <class R> class CPolyhedralDomain {
		typedef CPolyheder<R> Polyheder;
		//template <class G> friend bool operator == (const CPolyhedralDomain<G>& a, const CPolyhedralDomain<G>& b);
		//template <class G> friend bool operator == (const CPolyheder<G>& a, const CPolyheder<G>& b);
	private:
		int spaceDimension;
		//int parameterCount;
		vector<Polyheder* > polyhedra;

		void deletePolyhedra();
	protected:
	public:
		CPolyhedralDomain(int iSpaceDimension) : spaceDimension(iSpaceDimension) { };
		CPolyhedralDomain(const CPolyhedralDomain& iOriginal);

		void addPolyheder(const Polyheder& polyheder);
		const CPolyheder<R>& getPolyheder(int index) const { return *polyhedra[index]; }

		int getSpaceDimension() const { return spaceDimension; };
		int getPolyhederCount() const { return polyhedra.size(); };

		CPolyhedralDomain& operator = (const CPolyhedralDomain& other);

		CInterval<typename R::RationalType> asInterval() const;

		void projectFourierMotzkin(int projectionDimension);
		CPolyhedralDomain getProjectionFourierMotzkin(int projectionDimension) const;

		CFlat<CModule<R> > getAffineHull() const;
		CPolyheder<R> getConvexHull() const;
		CPolyheder<R> getMinimalPositiveHull() const;

		static CPolyhedralDomain<R> universe(int iSpaceDimension);

		//CPolyhedralDomain polyhedralLatticeDecomposition(const CLattice<R> core);

		operator string();

		bool isEmpty() const;
		bool isSpaceCompatibleWith(const CPolyhedralDomain<R>& other) const;

		string toStringHtml() const;
		string toString() const;
		void print() const;
		//int getParameterCount() const { return parameterCount; };

		CPolyhedralDomain<R> getNegative() const;

#ifdef PPL_LINK
		operator Parma_Polyhedra_Library::Powerset<ParmaCPolyhedron>() const;
#endif
	};

	template <class G> bool operator == (const CPolyhedralDomain<G>& a, const CPolyhedralDomain<G>& b);
	template <class G> bool operator != (const CPolyhedralDomain<G>& a, const CPolyhedralDomain<G>& b);
	template <class G> CPolyhedralDomain<G> operator && (const CPolyhedralDomain<G>& a, const CPolyhedralDomain<G>& b);
	template <class G> CPolyhedralDomain<G> operator && (const CPolyheder<G>& a, const CPolyhedralDomain<G>& b);
	template <class G> CPolyhedralDomain<G> operator && (const CPolyhedralDomain<G>& a, const CPolyheder<G>& b);
	template <class G> CPolyhedralDomain<G> operator || (const CPolyhedralDomain<G>& a, const CPolyhedralDomain<G>& b);

	template <class R>
	CPolyhedralDomain<R> CPolyhedralDomain<R>::universe(int iSpaceDimension) {
		CPolyhedralDomain<R> result = CPolyhedralDomain<R>(iSpaceDimension);

		result.addPolyheder(CPolyheder<R>::universe(iSpaceDimension));

		return result;
	}
}

#include "PolyhedralDomain.cpp"

#endif /* POLYHEDER_H_ */
