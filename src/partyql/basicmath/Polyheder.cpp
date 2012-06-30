#ifdef POLYHEDER_H_

#include "Polyheder.h"
#include "Lattice.h"
#include "Vektor.h"
#include "Flat.h"
#include "Cone.h"
#include "PolyhedralDomain.h"

#include "../../utils/Html.h"
#include "../core/util/Util.h"

#ifdef PPL_LINK
#include "ppl.h"
#endif

namespace AlgoTrans {
/*	template <class R> CPolyheder<R>::CPolyheder(int iSpaceDimension, vector<CVector<R>*> iTransformationVectors)
	: matrix(CPolyheder<R>(iSpaceDimension, iTransformationVectors))
	{
	}*/
	template <class R> CPolyheder<R>::CPolyheder(const CPolyheder<R>& iOriginal)
	{
		if (this != &iOriginal) *this = iOriginal;
	}

	template <class R> CPolyheder<R>& CPolyheder<R>::operator = (const CPolyheder<R>& other) {
		if (this != &other) matrix = other.getConstraintMatrix();

		return *this;
	}

	/*template <class G> CPolyheder<G> operator - (const CPolyheder<G>& a, const CPolyheder<G>& b) {
		//ASSERT(a.getSpaceDimension() == b.getSpaceDimension());

		CPolyheder<G> negB = CPolyheder<G>(b.getSpaceDimension());
		const CMatrix<G>& bConstr = b.getConstraintMatrix();
		for (int q = 0; q < b.getConstraintCount(); q++) {
			negB.addInequality((bConstr[q][0] + 1) << -bConstr[q].getSubVector(1, b.getSpaceDimension()));
		}

		return a && negB;
	}*/

	template <class R> CPolyheder<R>::operator CPolyhedralDomain<R>() const {
		CPolyhedralDomain<R> result = CPolyhedralDomain<R>(getSpaceDimension());

		result.addPolyheder(*this);

		return result;
	}

	template <class R> CPolyhedralDomain<R> CPolyheder<R>::getNegative() const {
		CPolyhedralDomain<R> result = CPolyhedralDomain<R>(getSpaceDimension());

		for (int q = 0; q < getConstraintCount(); q++) {
			result.addPolyheder(CMatrix<R>(-matrix[q] - CVector<R>::getUnitVector(getSpaceDimension() + 1, 0)));
		}

		return result;
	}

	template <class R> CPolyheder<R> operator + (const CPolyheder<R>& a, const CPolyheder<R>& b) {
		if (a.isEmpty()) return b;
		if (b.isEmpty()) return a;

		CCone<R> wA = CCone<R>(C, matrixToDescriptorSet(a.getFarkasPolyheder().getConstraintMatrix()));
		CCone<R> wB = CCone<R>(C, matrixToDescriptorSet(b.getFarkasPolyheder().getConstraintMatrix()));

		return CPolyheder<R>(DescriptorSetToMatrix((wA + wB).getDescriptors(C)));
	}

	template <class R>
	bool isNewlySuperfluous(CVector<R>& c, list<CVector<R>* > otherC) {
		typedef typename R::RationalType Rational;
		typedef typename list<CVector<R>* >::iterator listIt;

		CMatrix<R> (*U)(int dim) = &(CMatrix<R>::getUnitMatrix);
		CMatrix<R> (*Z)(int rows, int cols) = &(CMatrix<R>::getZeroMatrix);
		//CVector<R> (*UV)(int dim, int unitDim) = &(CVector<R>::getUnitVector);
		CVector<R> (*ZV)(int dim) = &(CVector<R>::getZeroVector);

		if (otherC.size() == 0) return ((c.getSubVector(1, c.getLength() - 1) == ZV(c.getLength() - 1)) && (c[0] >= 0));

		CIpProblem<R> ipp = CIpProblem<R>(otherC.size());

		CMatrix<R> consM = CMatrix<R>(c.getLength() - 1); // c == k . otherC ---- k >= 0
		consM.addRow(-c.getSubVector(1, c.getLength() - 1));
		for (listIt lit = otherC.begin(); lit != otherC.end(); ++lit) consM.addRow((**lit).getSubVector(1, c.getLength() - 1));

		CMatrix<R> domConM = consM.getTransposedMatrix() >> -consM.getTransposedMatrix() // c == k.otherC
		                     >> (Z(otherC.size(), 1) << U(otherC.size())); // k >= 0


		//domConM.print(); c.print();

		ipp.setDomainPolyheder(CFlat<CCone<R> >(1, CCone<R>::fromConstraints(domConM)));

		CVector<Rational>* solution = ipp.solve();

		bool nonSuperfluous = (solution != NULL);

		delete solution;

		return nonSuperfluous;
	}

	template <class R>
	bool isViolated(CVector<R>& c, const CPolyheder<R>& a) {
		typedef typename R::RationalType Rational;

		//CMatrix<R> (*U)(int dim) = &(CMatrix<R>::getUnitMatrix);
		CMatrix<R> (*Z)(int rows, int cols) = &(CMatrix<R>::getZeroMatrix);
		CVector<R> (*UV)(int dim, int unitDim) = &(CVector<R>::getUnitVector);
		//CVector<R> (*ZV)(int dim) = &(CVector<R>::getZeroVector);

		//if (otherC.size() == 0) return ((c.getSubVector(1, c.getLength() - 1) == ZV(c.getLength() - 1)) && (c[0] >= 0));

		CIpProblem<R> ipp = CIpProblem<R>(a.getSpaceDimension() + 1);

		CMatrix<R> aConstr = a.getConstraintMatrix();
		CMatrix<R> consM = CMatrix<R>(c.getLength() + 1);
		consM.addRow(-c.getSubVector(0, 1) << 1 << -c.getSubVector(1, c.getLength() - 1));
		consM.addRow(-(-c.getSubVector(0, 1) << 1 << -c.getSubVector(1, c.getLength() - 1)));
		consM >>= aConstr.getSubMatrixColumns(0, 1) << Z(aConstr.getRowCount(), 1)
					<< aConstr.getSubMatrixColumns(1, aConstr.getColumnCount() - 1);
		consM.addRow((c.getSubVector(0, 1) + UV(1, 0)) << 0 << c.getSubVector(1, c.getLength() - 1));
		/* k = c.x
		 * A.x >= 0
		 * c.x >= -1
		 */

		//domConM.print(); c.print();

		ipp.setDomainPolyheder(CFlat<CCone<R> >(1, CCone<R>::fromConstraints(consM)));

		CVector<Rational>* solution = ipp.solve();

		bool isViolated = (solution != NULL) && ((*solution)[0] < CRational<R>(0));

		delete solution;

		return isViolated;
	}

	template <class R>
	bool isNonMinimal(CVector<R>& c, list<CVector<R>* > otherC) {
		typedef typename R::RationalType Rational;
		typedef typename list<CVector<R>* >::iterator listIt;

		CMatrix<R> (*U)(int dim) = &(CMatrix<R>::getUnitMatrix);
		CMatrix<R> (*Z)(int rows, int cols) = &(CMatrix<R>::getZeroMatrix);
		CVector<R> (*UV)(int dim, int unitDim) = &(CVector<R>::getUnitVector);
		CVector<R> (*ZV)(int dim) = &(CVector<R>::getZeroVector);

		if (otherC.size() == 0) return ((c.getSubVector(1, c.getLength() - 1) == ZV(c.getLength() - 1)) && (c[0] >= 0));

		CMatrix<R> consM = CMatrix<R>(c.getLength()); // c == k . otherC ---- k >= 0
		consM.addRow(-c);
		for (listIt lit = otherC.begin(); lit != otherC.end(); ++lit) consM.addRow(**lit);
		consM.addRow(UV(c.getLength(), 0));

		CMatrix<R> domConM = consM.getTransposedMatrix() >> -consM.getTransposedMatrix() // c == k.otherC + k0
		                     >> (Z(otherC.size(), 1) << U(otherC.size()) << Z(otherC.size(), 1)); // k >= 0 //&& k0 >= 0


		//domConM.print(); c.print();

		ASSERT(false);
		/* XXX: readd ParmaMIPProblem plpp = ParmaMIPProblem((ParmaConstraintSystem) domConM,
				                (ParmaLinearExpression) UV(otherC.size() + 1, 1),
				                Parma_Polyhedra_Library::MINIMIZATION);*/

//		plpp.solve();

		//print();

		return false;
		// XXX: return plpp.is_satisfiable();
	}

	template <class R> CPolyheder<R> CPolyheder<R>::getAMinimalPositiveHull() const {
		typedef typename list<CVector<R>* >::iterator listIt;
		list<CVector<R>* > consList = list<CVector<R>* >();
		for (int q = 0; q < matrix.getRowCount(); q++) consList.push_back(&matrix[q]);

		for (listIt lit = consList.begin(); lit != consList.end(); ) {
			list<CVector<R>* > othC;

			for (listIt slit = consList.begin(); slit != consList.end(); ++slit) {
				if (slit != lit) othC.push_back(*slit);
			}
			if (isNonMinimal<R>(**lit, othC)) {
				listIt pit = lit;
				++lit;
				consList.erase(pit);
			} else ++lit;
		}

		CMatrix<R> resultCons = CMatrix<R>(1 + getSpaceDimension());
		typedef typename list<CVector<R>* >::iterator listIt;
		for (listIt lit = consList.begin(); lit != consList.end(); ++lit) resultCons.addRow(**lit);

		return (resultCons.getRowCount() == 0) ? CPolyheder<R>::universe(getSpaceDimension()) : CPolyheder<R>(resultCons);
	}

	template <class G>
	CPolyheder<G> getCompatiblePolyhedron(const CPolyheder<G>& a, const CPolyheder<G>& b) {
		CMatrix<G> cons = a.getConstraintMatrix();

		CMatrix<G> nonViolated = CMatrix<G>(cons.getColumnCount());
		for (int q = 0; q < cons.getRowCount(); q++) {
			if (!isViolated(cons[q], b)) nonViolated >>= cons[q];
		}

		return CPolyheder<G>(nonViolated);
	}

	template <class G> CPolyheder<G> operator || (const CPolyheder<G>& a, const CPolyheder<G>& b) {
		if (a.isEmpty()) return b;
		if (b.isEmpty()) return a;

		CPolyheder<G> aPrep = getCompatiblePolyhedron(a, b).getAMinimalPositiveHull();
		aPrep.print();
		CPolyheder<G> bPrep = getCompatiblePolyhedron(b, aPrep).getAMinimalPositiveHull();
		bPrep.print();

		/*typedef typename list<CVector<G>* >::iterator listIt;
		CMatrix<G> cons = aPrep.getConstraintMatrix() >> bPrep.getConstraintMatrix();

		list<CVector<G>* > consList = list<CVector<G>* >();
		for (int q = 0; q < cons.getRowCount(); q++) consList.push_back(&cons[q]);

		for (listIt lit = consList.begin(); lit != consList.end(); ) {
			list<CVector<G>* > othC;

			for (listIt slit = consList.begin(); slit != consList.end(); ++slit) {
				if (slit != lit) othC.push_back(*slit);
			}
			if (isNewlySuperfluous<G>(**lit, othC)) {
				listIt pit = lit;
				++lit;
				consList.erase(pit);
			} else ++lit;
		}

		CMatrix<G> resultCons = CMatrix<G>(1 + a.getSpaceDimension());
		typedef typename list<CVector<G>* >::iterator listIt;
		for (listIt lit = consList.begin(); lit != consList.end(); ++lit) resultCons.addRow(**lit);*/

		return CPolyheder<G>(aPrep.getConstraintMatrix() >> bPrep.getConstraintMatrix()).getAMinimalPositiveHull();
	}

	template <class G> CPolyheder<G> operator && (const CPolyheder<G>& a, const CPolyheder<G>& b) {
		if (a.getConstraintMatrix().getRowCount() == 0) return b;
		if (b.getConstraintMatrix().getRowCount() == 0) return a;

		return CPolyheder<G>(a.getConstraintMatrix() >> b.getConstraintMatrix()).getSimplifiedPolyheder();
	}

	template <class G> CPolyhedralDomain<G> operator - (const CPolyheder<G>& a, const CPolyheder<G>& b) { // (Reflexive) SubSetOf
		return ((CPolyhedralDomain<G>) a) && b.getNegative();
	}

	template <class G> bool operator <= (const CPolyheder<G>& a, const CPolyheder<G>& b) { // (Reflexive) SubSetOf
		return (a - b).isEmpty();
	}

	template <class G> bool operator == (const CPolyheder<G>& a, const CPolyheder<G>& b) {
		return (a <= b) && (b <= a);
	}

	template <class G> bool operator != (const CPolyheder<G>& a, const CPolyheder<G>& b) {
		return !(a == b);
	}

	template <class R> bool CPolyheder<R>::isEmpty() const {
		ASSERT(false);
		/* XXX:
		ParmaMIPProblem plpp = ParmaMIPProblem(getSpaceDimension(), (ParmaConstraintSystem) getConstraintMatrix(),
				                (ParmaLinearExpression) CVector<R>::getUnitVector(getSpaceDimension() + 1, 1),
				                Parma_Polyhedra_Library::MINIMIZATION);

		return !plpp.is_satisfiable();*/

		return false;

		/* Correct, but slower (using projection):
		CCone<R> p = CCone<R>::fromConstraints(getFarkasPolyheder().getConstraintMatrix());

		for (int q = 0; q < getSpaceDimension(); q++) {//p.print();
			p = p.getProjectionFourierMotzkin(1);
		}//p.print();

		CMatrix<R> m = p.getConstraintMatrix();

		for (int r = 0; r < m.getRowCount(); r++) {
			if (m[r][0] < 0) return true;
		}

		return false;*/
	}

	template <class R> CInterval<typename R::RationalType> CPolyheder<R>::asInterval() const {
		typedef typename R::RationalType Rat;

		ASSERT(getSpaceDimension() == 1);

		bool lowerBounded = false, upperBounded = false;
		Rat lowerBound, upperBound;
		const CMatrix<R>& X = this->getConstraintMatrix();
		for (int r = 0; r < getConstraintCount(); r++) {
			if (X(r, 1) != 0) {
				Rat newCandidate = -X(r, 0)/X(r, 1);
				if (X(r, 1) > 0) {
					lowerBound = !lowerBounded ? newCandidate : ((newCandidate > lowerBound) ? newCandidate : lowerBound);
					lowerBounded = true;
				} else {
					upperBound = !upperBounded ? newCandidate : ((newCandidate < upperBound) ? newCandidate : upperBound);
					upperBounded = true;
				}
			}
		}

		return CInterval<Rat>(lowerBounded, lowerBound, upperBounded, upperBound);
	}

	template <class R> CInterval<typename R::RationalType> CPolyheder<R>::asIntervalAlongDimension(int dim) const {
		typedef typename R::RationalType Rat;

		CPolyheder<R> imdPoly = *this;
		for (int q = 0; q < dim; q++) {
			imdPoly = imdPoly.getProjectionFourierMotzkin(0);
		}
		for (int q = imdPoly.getSpaceDimension() - 1; q > 0; q--) {
			imdPoly = imdPoly.getProjectionFourierMotzkin(q);
		}

		return imdPoly.asInterval();
	}

	template <class R> CPolyheder<R> CPolyheder<R>::getProjectionFourierMotzkin(int projectionDim) const {
		ASSERT((projectionDim >= 0) && (projectionDim < getSpaceDimension())); // We're allowing projection of the constant too for cone-projections...

		// Use cone's projection
		CCone<R> cone = CFlat<CCone<R> >(1,
				CCone<R>(C, matrixToDescriptorSet(getFarkasPolyheder().getConstraintMatrix()))).getProjectionFourierMotzkin(projectionDim + 1);

		return CPolyheder<R>(cone.matrix());
	}

	template <class R> CPolyheder<R> CPolyheder<R>::getProjectionJones(int projectionDim) const {
		ASSERT((projectionDim >= 0) && (projectionDim < getSpaceDimension())); // We're allowing projection of the constant too for cone-projections...

		// Use cone's projection
		CCone<R> cone = CFlat<CCone<R> >(1,
				CCone<R>(C, matrixToDescriptorSet(getFarkasPolyheder().getConstraintMatrix()))).getProjectionJones_ThroughPIP(getInterval(projectionDim + 1, 1));

		return CPolyheder<R>(cone.matrix());
	}

	template <class R> CPolyheder<R> CPolyheder<R>::getFarkasPolyheder() const {
		return getConstraintMatrix()
		       >> CVector<R>::getUnitVector(getSpaceDimension() + 1, 0);
	}

	/*template <class R> bool CPolyheder<R>::isEmpty() const {
		CCone<R> homogCone = CCone<R>::fromConstraints(getConstraintMatrix()
				                                       >> CVector<R>::getUnitVector(getSpaceDimension() + 1, 0));

		return CFlat<CCone<R> >(1, homogCone).isPoint();
	}*/

	template <class R> CFlat<CModule<R> > CPolyheder<R>::getAffineHull_LimsLinearisation() const {
		// Essentially Lim's Linearisation method (without equality optimizations)
		CMatrix<R> (*U)(int dim) = &(CMatrix<R>::getUnitMatrix);
		CMatrix<R> (*Z)(int rows, int cols) = &(CMatrix<R>::getZeroMatrix);
		//CVector<R> (*UV)(int dim, int unitDim) = &(CVector<R>::getUnitVector);

		CMatrix<R> farkasT = getFarkasPolyheder().getConstraintMatrix().getTransposedMatrix();
		int farCount = farkasT.getColumnCount();

		CCone<R> p = CCone<R>(C, matrixToDescriptorSet( (farkasT << -U(getSpaceDimension() + 1))
		                  >> -(farkasT << -U(getSpaceDimension() + 1))
		                  >>  (U(farCount)                  << Z(farCount, getSpaceDimension() + 1))));

		for (int q = 0; q < farCount; q++) {
			p = p.getProjectionFourierMotzkin(0, (q == 0));
		}
		CMatrix<R> mCon = !DescriptorSetToMatrix(p.getDescriptors(C));

		return CFlat<CModule<R> >(1, CModule<R>(C, matrixToDescriptorSet(mCon)));
	}

	template <class R> CFlat<CModule<R> > CPolyheder<R>::getAffineHull() const {
		CMatrix<R> (*U)(int dim) = &(CMatrix<R>::getUnitMatrix);
		CMatrix<R> (*Z)(int rows, int cols) = &(CMatrix<R>::getZeroMatrix);
		CVector<R> (*CV)(int dim, R constant) = &(CVector<R>::getConstantVector);

		if (isEmpty()) return CFlat<CModule<R> >(1, CModule<R>(getSpaceDimension() + 1, G));

		CMatrix<R> resultCons = CMatrix<R>(getSpaceDimension() + 1);
		CMatrix<R> consM = getFarkasPolyheder().getConstraintMatrix();
		//consM.print();
		for (int s = getSpaceDimension(); s > 0; s--) {
			CMatrix<R> consT = consM.getTransposedMatrix();
			CMatrix<R> domCons =
				(Z(consT.getRowCount(), 1) << consT << consT)
				>> -(Z(consT.getRowCount(), 1) << consT << consT)
				>> (Z(2*consT.getColumnCount(), 1) << U(2*consT.getColumnCount()))
				>> (CV(1, -1) << CV(consT.getColumnCount(), 1) << CV(consT.getColumnCount(), 0))
				>> (CV(1, -1) << CV(consT.getColumnCount(), 0) << CV(consT.getColumnCount(), 1));
			//domCons.print();

			ASSERT(false);
			/* XXX:
			ParmaMIPProblem plpp = ParmaMIPProblem(domCons.getColumnCount() - 1, (ParmaConstraintSystem) domCons,
					                (ParmaLinearExpression) CVector<R>::getUnitVector(domCons.getColumnCount(), 1),
					                Parma_Polyhedra_Library::MINIMIZATION);
			plpp.solve();
			if (!plpp.is_satisfiable()) break;


			CVector<R> sol = CVector<R>::fromParmaGenerator(plpp.optimizing_point());
			ASSERT(sol[0] == (R) 1);
			sol = sol.getSubVector(1, sol.getLength() - 1);
			//sol.print();
			CVector<R> kerV = CV(consM.getColumnCount(), 0);
			for (int q = 0; q < consM.getRowCount(); q++) kerV += sol[q] * consM[q];
			//kerV.print();
			resultCons.addRow(kerV);

			// Reduce remaining constraint matrix
			if (s > 1) { // unless we do not need it anymore
				R normSqKerV = kerV.normSquare();
				//normSqKerV.print();
				CMatrix<R> consMNew = CMatrix<R>(consM.getColumnCount());
				for (int q = 0; q < consM.getRowCount(); q++) {
					//consM[q].print();
					R dotProd = consM[q].innerProduct(kerV);
					//dotProd.print();
					if (dotProd != 0) {
						R g = gcd(dotProd, normSqKerV);
						//g.print();
						R mK = dotProd / g;
						//mK.print();
						R mC = normSqKerV / g;
						//mC.print();
						CVector<R> newRow = mC * consM[q] - mK * kerV;
						if (!(newRow.isZero())) consMNew.addRow(newRow);
						//newRow.print();
					} else consMNew.addRow(consM[q]);
				}
				consM = consMNew;
			}
			*/
		}

		//resultCons.print();

		return CFlat<CModule<R> >(1, CModule<R>(C, matrixToDescriptorSet(resultCons)));
	}

	template <class R> CPolyheder<R> CPolyheder<R>::polyhedralLatticeDecomposition(const CLattice<R> core) {
		const CMatrix<R>& X = this->getConstraintMatrix();

		CMatrix<R> coreMatrix = CMatrix<R>(core.getMatrix());
		CHermiteNormalFormReducer<R> reducer = CHermiteNormalFormReducer<R>(coreMatrix);
		reducer.reduce();

		vector<int> usefulPivots, spaceDimensions;
		vector<R> multipliers;
		for (unsigned int p = 0; p < reducer.pivotColumns.size(); p++) {
			if (coreMatrix[p][reducer.pivotColumns[p]] != 1) {
				usefulPivots.push_back(reducer.pivotColumns[p]);
				multipliers.push_back(coreMatrix[p][reducer.pivotColumns[p]]);
			}
		}

		int newSpaceDim = getSpaceDimension() + multipliers.size();
		CPolyheder<R> newPolyh = CPolyheder<R>(newSpaceDim);
		for (int r = 0; r < getConstraintCount(); r++) {
			unsigned int nextPivot = 0;
			CVector<R> row;
			for (int c = 0; c < matrix.getColumnCount(); c++) {
				int nextPivIx = (nextPivot < usefulPivots.size()) ? usefulPivots[nextPivot] : -1;
				if (nextPivIx == c) {
					row.appendElement(multipliers[nextPivIx] * X(r, c));
					row.appendElement(X(r, c));
					nextPivot++;
				} else {
					row.appendElement(X(r, c));
				}
			}
			newPolyh.addInequality(row);
		}

		// Add elementary parallelotope constraints
		int nextPivot = 0;
		int nC = 0;
		for (int c = 0; c < newSpaceDim; c++) {
			int nextPivIx = (nextPivot < (int) usefulPivots.size()) ? usefulPivots[nextPivot] : -1;
			if (nextPivIx == c) {
				nC++;
				newPolyh.addInequality(CVector<R>::getUnitVector(newSpaceDim + 1, nC));
				newPolyh.addInequality((multipliers[nextPivIx] - 1) * CVector<R>::getUnitVector(newSpaceDim + 1, newSpaceDim) - CVector<R>::getUnitVector(newSpaceDim + 1, nC));
				nC++;
				nextPivot++;
			} else nC++;
		}

		return newPolyh;
	}

	template <class R> CPolyheder<R> CPolyheder<R>::getDual() const {
		CCone<R> w = CCone<R>(C, matrixToDescriptorSet(getFarkasPolyheder().getConstraintMatrix()));

		return CPolyheder<R>(DescriptorSetToMatrix((!w).getDescriptors(C)));
	}

	template <class R> CPolyheder<R> CPolyheder<R>::fromConstraints(const CMatrix<R>& equalities, const CMatrix<R>& inequalities) {
		return CPolyheder(equalities >> -equalities >> inequalities);
	}

	template <class R> CPolyheder<R> CPolyheder<R>::fromInequalities(const CMatrix<R>& inequalities) {
		return CPolyheder(inequalities);
	}

	template <class R> CPolyheder<R> CPolyheder<R>::fromEqualities(const CMatrix<R>& equalities) {
		return fromInequalities(equalities >> -equalities);
	}

	template <class R> CPolyheder<R> CPolyheder<R>::getSimplifiedPolyheder() const {
		return CPolyheder<R>::fromInequalities(CCone<R>(C, matrixToDescriptorSet(getFarkasPolyheder().getConstraintMatrix())).getSimplifiedConeFast().matrix());
	}

#ifdef PPL_LINK
	template <class R>
	CPolyheder<R>::operator ParmaCPolyhedron() const {
		return ParmaCPolyhedron((ParmaConstraintSystem) matrix);
	}

	template <class R>
	CPolyheder<R> CPolyheder<R>::fromParmaCPolyhedron(const Parma_Polyhedra_Library::C_Polyhedron& polyhedron) {
		return CPolyheder<R>(CMatrix<R>::fromParmaConstraintSystem(polyhedron.constraints()));
	}
#endif

	template <class R>
	string CPolyheder<R>::toStringHtml() const {
		return Html::wrapper("Polyheder", matrix.toStringHtml());
	}
}

#endif /* POLYHEDER_CPP_ */
