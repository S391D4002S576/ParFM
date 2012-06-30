#ifndef POLYHEDER_H_
#define POLYHEDER_H_

#include <vector>
using std::vector;

#include <list>
using std::list;

#include <algorithm>
using std::min;

#include "Vektor.h"
#include "Matrix.h"
#include "Module.h"
#include "scalar/Interval.h"

#include <string>
using std::string;

#include "../cute/cute.h"

//#include <ppl.hh>

namespace AlgoTrans {
#ifdef PPL_LINK
	typedef class Parma_Polyhedra_Library::C_Polyhedron ParmaCPolyhedron;
#endif

	template <class R> class CLattice;
	template <class R> class CPolyhedralDomain;
	template <class R> class CFlat;

	template <class R> class CPolyheder {
		template <class G> friend CPolyheder<G> operator || (const CPolyheder<G>& a, const CPolyheder<G>& b);
	private:
		CMatrix<R> matrix;
	protected:

	public:
		CPolyheder() { }
		CPolyheder(int iSpaceDimension) : matrix(CMatrix<R>(iSpaceDimension + 1)) { }
		CPolyheder(const CMatrix<R>& iMatrix) : matrix(iMatrix) { }
		CPolyheder(const CPolyheder& iOriginal);

		const CMatrix<R>& getConstraintMatrix() const { return matrix; }

		void addInequality(const CVector<R>& inequalityVector) { matrix.addRow(inequalityVector); }
		void addEquality(const CVector<R>& equalityVector) { addInequality(-equalityVector); addInequality(equalityVector); }

		static CPolyheder<R> fromConstraints(const CMatrix<R>& equalities, const CMatrix<R>& inequalities);
		static CPolyheder<R> fromInequalities(const CMatrix<R>& inequalities);
		static CPolyheder<R> fromEqualities(const CMatrix<R>& equalities);

		static CPolyheder<R> source(int iSpaceDimension) { return CPolyheder<R>::fromInequalities(-CVector<R>::getUnitVector(iSpaceDimension + 1, 0)); };
		static CPolyheder<R> universe(int iSpaceDimension) { return CPolyheder<R>::fromInequalities(CVector<R>::getUnitVector(iSpaceDimension + 1, 0)); };

		int getSpaceDimension() const { return matrix.getColumnCount() - 1; };
		int getConstraintCount() const { return matrix.getRowCount(); };

		CPolyheder<R> decompose(vector<R> coreLatticeMultipliers, vector<int> decomposedDims);

		CPolyheder<R>& operator = (const CPolyheder<R>& other);

		CInterval<typename R::RationalType> asInterval() const;
		CInterval<typename R::RationalType> asIntervalAlongDimension(int dimension) const;

		CPolyheder<R> getProjectionFourierMotzkin(int projectionDimension) const;
		CPolyheder<R> getProjectionJones(int projectionDimension) const;

		CPolyheder<R> getFarkasPolyheder() const;

		CFlat<CModule<R> > getAffineHull_LimsLinearisation() const;
		CFlat<CModule<R> > getAffineHull() const;
		bool isEmpty() const;

		CPolyheder<R> polyhedralLatticeDecomposition(const CLattice<R> core);

		string toStringHtml() const;
		string toString() const { return "Polyheder(" + matrix.toString() + ")"; }
		void print() const { printf((toString() + "\n").c_str()); }

		CPolyhedralDomain<R> getNegative() const;
		CPolyheder<R> getDual() const;

		CPolyheder<R> getSimplifiedPolyheder() const;
		CPolyheder<R> getAMinimalPositiveHull() const;

		operator CPolyhedralDomain<R>() const;

#ifdef PPL_LINK
		operator Parma_Polyhedra_Library::C_Polyhedron() const;
		static CPolyheder<R> fromParmaCPolyhedron(const Parma_Polyhedra_Library::C_Polyhedron& polyhedron);
#endif
	};

	template <class G> bool operator == (const CPolyheder<G>& a, const CPolyheder<G>& b);
	template <class G> CPolyheder<G> operator && (const CPolyheder<G>& a, const CPolyheder<G>& b);
	template <class G> CPolyheder<G> operator + (const CPolyheder<G>& a, const CPolyheder<G>& b); // Convex hull
	template <class G> CPolyheder<G> operator || (const CPolyheder<G>& a, const CPolyheder<G>& b); // Convex hull
}

#include "Polyheder.cpp"

#endif /* POLYHEDER_H_ */
