#ifndef CONE_H_
#define CONE_H_

#include <vector>
using std::vector;

#include <list>
using std::list;

#include <algorithm>
using std::min;

#include "Matrix.h"
#include "Hadron.h"
#include "VectorDescriptor.h"

#include "IntegerProgramming.h"
//#include "VectorDescriptor.h"

#include "../core/util/Time.h"
#include "../core/util/Util.h"

namespace AlgoTrans {
	template <class R> struct VectorRedundanceStruct {
		CVector<R>* v;
		bool certainlyNonRedundant;

		VectorRedundanceStruct(CVector<R>* iV, bool iCertainlyNonRedundant = false)
		: v(iV), certainlyNonRedundant(iCertainlyNonRedundant) { }
	};

	template <class R> class CCone : public CHadron<CUnidirectionalVectorDescriptorSet<R> > {
		typedef CUnidirectionalVectorDescriptorSet<R> DescriptorSet;

		template <class G> friend bool operator == (const CCone<G>& a, const CCone<G>& b);

		list<CVector<R>* > retainNonRedundant(list<CVector<R>* > otherC) const;
		bool isRedundant(CVector<R>& c, list<CVector<R>* > otherC) const;
		list<VectorRedundanceStruct<R>* > retainNonRedundant(list<VectorRedundanceStruct<R>* > otherC) const;
	protected:

	public:
		const CMatrix<R> matrix();

		CCone(int iSpaceDimension, HadronDescription description) : CHadron<DescriptorSet>(iSpaceDimension, description) { };
		CCone(const CCone<R>& iOriginal) : CHadron<DescriptorSet>(iOriginal) { };
		CCone(HadronDescription description, const DescriptorSet& descriptors) : CHadron<DescriptorSet>(description, descriptors) { };
		CCone(const CHadron<DescriptorSet>& iOriginal) : CHadron<DescriptorSet>(iOriginal) { };

		CCone<R> getProjectionFourierMotzkin(int projectionDim, bool preSimplify = true);
		CCone<R> getProjectionJones_ThroughPIP(std::vector<int> removedDims);

		CCone<R> getReflexiveTransitiveClosure(int spaceDim);

		CCone<R> getSimplifiedConeFast();
		CCone<R> getSimplifiedConeFast_Ex(int numberOfCertainlyNonRedundant);

		CCone<R> getDualizedCone();

#ifdef PPL_LINK
		CCone<R> getProjectionThroughVertexEnumeration(std::vector<int> removedDims);
		CCone<R> getDualizedCone_Parma();
#endif

		//static std::vector<int> getDimSetComplement(std::vector<int> dims, int dimension);
		//static std::vector<int> getInterval(int start, int length);
		//string toStringHtml() const { return ("Cone(" + constraintMatrix.toStringHtml() + ")"); };

		//CCone<R>& operator = (const CCone<R>& other);
	};

	/*template <class R> CCone<R> operator + (CCone<R> a, CCone<R> b) {
		return CCone<R>(((CHadron<CUnidirectionalVectorDescriptorSet<R> >) a) + ((CHadron<CUnidirectionalVectorDescriptorSet<R> >) b));
	}*/

	template <class R> CCone<R> operator * (CCone<R> a, CCone<R> b) {
		return CCone<R>(((CHadron<CUnidirectionalVectorDescriptorSet<R> >) a) * ((CHadron<CUnidirectionalVectorDescriptorSet<R> >) b));
	}

	template <class R> CCone<R> meet(CCone<R>& a, CCone<R>& b) {
		return CCone<R>(meet((CHadron<CUnidirectionalVectorDescriptorSet<R> >&) a), ((CHadron<CUnidirectionalVectorDescriptorSet<R> >&) b));
	}

	template <class R> CCone<R> join(CCone<R>& a, CCone<R>& b) {
		return CCone<R>(join((CHadron<CUnidirectionalVectorDescriptorSet<R> >&) a), ((CHadron<CUnidirectionalVectorDescriptorSet<R> >&) b));
	}

	/*template <class R> CCone<R> operator !(CCone<R>& a) {
		return CCone<R>(!((CHadron<CUnidirectionalVectorDescriptorSet<R> >&) a));
	}*/
#ifdef PPL_LINK
	typedef class Parma_Polyhedra_Library::MIP_Problem ParmaMIPProblem;
	typedef class Parma_Polyhedra_Library::Constraint_System ParmaConstraintSystem;
#endif

	static double timeConvexHull_J, timeConvexHull_FM /*, timeConvexHull_VE , timeConvexHull_Dual*/;

	template <class R> CCone<R> operator + (CCone<R> a, CCone<R> b) { ASSERT(a.getDimension() == b.getDimension());
		double tPrev;

		int d = a.getDimension();

		CMatrix<R> (*U)(int dim) = &(CMatrix<R>::getUnitMatrix);
		CMatrix<R> (*Z)(int rows, int cols) = &(CMatrix<R>::getZeroMatrix);

		CMatrix<R> A = DescriptorSetToMatrix(a.getDescriptors(C));
		CMatrix<R> B = DescriptorSetToMatrix(b.getDescriptors(C));

		//A.print();
		//B.print();

		CCone<R> res2 = CCone<R>(C, matrixToDescriptorSet(
			(A << Z(A.getRowCount(), 2*d))
			>> (Z(B.getRowCount(), d) << B << Z(B.getRowCount(), d))
			>> (U(d) << U(d) << -U(d))
			>> -(U(d) << U(d) << -U(d))
		));
		CCone<R> res4 = CCone<R>(C, matrixToDescriptorSet(
			(A << Z(A.getRowCount(), d))
			>> ((-B) << B)
		));

		/*CMatrix<R> G = A.getTransposedMatrix();
		CMatrix<R> H = B.getTransposedMatrix();
		CCone<R> res3 = CCone<R>(C, matrixToDescriptorSet(
			(Z(d, d) << -U(d) << G << Z(d, H.getColumnCount()))
			>> (Z(d, d) << -U(d) << Z(d, G.getColumnCount()) << H)
			>> (Z(G.getColumnCount(), d) << Z(G.getColumnCount(), d) << U(G.getColumnCount()) << Z(G.getColumnCount(), H.getColumnCount()))
			>> (Z(H.getColumnCount(), d) << Z(H.getColumnCount(), d) << Z(H.getColumnCount(), G.getColumnCount()) << U(H.getColumnCount()))
			>> (U(d) << Z(d, G.getColumnCount() + H.getColumnCount()) << -U(d))
		));*/
		//res2.print();

		tPrev = rtclock();
		CCone<R> res3 = res4;
		for (int q = d - 1; q >= 0; q--) {
			//DescriptorSetToMatrix(res2.getDescriptors(C)).print();
			//printf("%d:", res3.getDescriptors(C).getSize());
			res3 = res3.getProjectionFourierMotzkin(q, q == d - 1);
		}
		//printf("\n");
		double timeConvexHull_FMDuration = (rtclock() - tPrev);
		timeConvexHull_FM += (rtclock() - tPrev);
		//res3.print();

#ifdef PPL_LINK
		tPrev = rtclock();
		CCone<R> res6 = res2.getProjectionThroughVertexEnumeration(getInterval(0, 2*d));
		double timeConvexHull_VEDuration = (rtclock() - tPrev);
		timeConvexHull_VE += (rtclock() - tPrev);
#endif
		//res3.print();
		//res6.print();
		//ASSERT(res6 == res3);

		tPrev = rtclock();
		CCone<R> res5 = res2.getProjectionJones_ThroughPIP(getInterval(0, 2*d));
		//ASSERT(res5 == res3);

		/*CVector<R> (*UV)(int dim, int unitDim) = &(CVector<R>::getUnitVector);
		CVector<R> (*ZV)(int dim) = &(CVector<R>::getZeroVector);

		CMatrix<R> res2udM = (A << Z(A.getRowCount(), 2*d))
							 >> (Z(B.getRowCount(), d) << B << Z(B.getRowCount(), d));
		CMatrix<R> res2bdM = (U(d) << U(d) << -U(d));

		// calc [C b]a == c a
		CMatrix<R> CM2 = CMatrix<R>(2 + res2udM.getColumnCount());
		for (int q = 0; q < res2udM.getRowCount(); q++) {
			CM2 = CM2 >> (UV(2, 1) << res2udM[q]);
			CM2 = CM2 >> (ZV(2) << res2udM[q]);
		}
		CM2 = CM2 >> -(UV(CM2.getColumnCount(), 1));

		//CM2.print();
		CPipProblem<R> pp = CPipProblem<R>(1 + 2*d, d);
		pp.integerSolution = false;
		pp.simplify = false;
		pp.domConstrFlat = CFlat<CModule<R> >(1, CModule<R>(C, matrixToDescriptorSet(Z(res2bdM.getRowCount(), 2) << res2bdM)));
		pp.domConstrPolyheder = CFlat<CCone<R> >(1, CCone<R>(C, matrixToDescriptorSet(CM2)));

		CPipQuast<R>* solution = pp.solve();

		ASSERT(solution != NULL);
		//if (solution != NULL) solution->print();

		CMatrix<R> newCs = extractConstraintMatrixFromPipQuast(d, solution);
		CCone<R> res5 = CCone<R>(C, matrixToDescriptorSet(newCs)).getSimplifiedConeFast();*/

		double timeConvexHull_JDuration = (rtclock() - tPrev);
		timeConvexHull_J += (rtclock() - tPrev);

		/*tPrev = rtclock();
		CCone<R> res6 = (CCone<R>) ((CHadron<CUnidirectionalVectorDescriptorSet<R> >) a
		                 + (CHadron<CUnidirectionalVectorDescriptorSet<R> >) b);
		timeConvexHull_Dual += (rtclock() - tPrev);*/

		printf("FM  = %0.6lfs [%0.6lfs]\n", timeConvexHull_FM, timeConvexHull_FMDuration);
		printf("J = %0.6lfs [%0.6lfs]\n", timeConvexHull_J, timeConvexHull_JDuration);
#ifdef PPL_LINK
		printf("VE = %0.6lfs [%0.6lfs]\n", timeConvexHull_VE, timeConvexHull_VEDuration);
#endif
		//printf("Dual = %0.6lfs\n", timeConvexHull_Dual);


		/*CCone<R> res1 = CCone<R>(G, a.getDescriptors(G) + b.getDescriptors(G));
		res1.print();
		res2.print();*/

		return res3;
	}

	template <class G> bool operator == (const CCone<G>& a, const CCone<G>& b) {
		ASSERTM("Not yet implemented", false);

		return false;
	}

	template <class R> CVector<R> _RVectorFromRationalVector(CVector<typename R::RationalType>& rationalVector) {
		R denLCM = R::getOne();

		for (int q = 0; q < rationalVector.getLength(); q++) {
			denLCM = lcm(denLCM, rationalVector[q].denominator);
		}

		CVector<R> result = CVector<R>();
		for (int q = 0; q < rationalVector.getLength(); q++) {
			result.appendElement(rationalVector[q].numerator * (denLCM / rationalVector[q].denominator));
		}

		return result;
	}

	template <class R> CMatrix<R> extractConstraintMatrixFromPipQuast(int dim, CPipQuast<R>* quast, std::string indent = "") {
		if (quast != NULL) {
			if (quast->condition != NULL) {
				//std::cout << indent + "Condition: " + _RVectorFromRationalVector<R>(*(quast->condition)).toString() + "\n";
				CMatrix<R> nonNeg = extractConstraintMatrixFromPipQuast(dim, quast->nonNegativeCondQuast, indent + ": ");
				//std::cout << indent + ")(\n";
				CMatrix<R> neg = extractConstraintMatrixFromPipQuast(dim, quast->negativeCondQuast, indent + ": ");
				//std::cout << indent + ")\n";
				return nonNeg >> neg;
			} else {
				//std::cout << indent + "Solution:";
				CMatrix<R> result = CMatrix<R>(dim);
				if (quast->solution != NULL) {
					//(_RVectorFromRationalVector<R>((*quast->solution)[0])).print();
					result = result >> -(_RVectorFromRationalVector<R>((*quast->solution)[0]).getSubVector(0, dim));
				} //else printf("[NULL]\n");
				return result;
			}
		} else {
			//std::cout << indent + "[Empty]";
			return CMatrix<R>(dim);
		}
	}

#ifdef PPL_LINK
	template <class R> CCone<R> CCone<R>::getProjectionThroughVertexEnumeration(std::vector<int> removedDims) {
		std::vector<int> retainedDims = getDimSetComplement(removedDims, this->getDimension());
		CMatrix<R> (*U)(int dim) = &(CMatrix<R>::getUnitMatrix);
		//CVector<R> (*UV)(int dim, int unitDim) = &(CVector<R>::getUnitVector);
		//CVector<R> (*CV)(int dim, R constant) = &(CVector<R>::getConstantVector);
		//CVector<R> (*ZV)(int dim) = &(CVector<R>::getZeroVector);

		CMatrix<R> Y = matrix();
		CMatrix<R> d = CMatrix<R>(removedDims.size()); // Fill a matrix with the constraint coefficients corresponding to removed dimensions
		ITT(std::vector<int>, i, removedDims) d = d << Y.getColumn(*i);

		CMatrix<R> c = CMatrix<R>(retainedDims.size()); // Fill a matrix with the constraint coefficients corresponding to retained dimensions
		ITT(std::vector<int>, i, retainedDims) c = c << Y.getColumn(*i);

		//c.print();
		//d.print();

		CMatrix<R> cDual = !(c.getTransposedMatrix());
		//cDual.print();

		d = d.getTransposedMatrix();
		d.reduceToHermiteNormalForm();

		CMatrix<R> Z = d
		               >> -d
		               //>> cDual
		               //>> -cDual
		               >> U(d.getColumnCount());
		CCone<R> ZCone = CCone<R>(C, matrixToDescriptorSet(Z));
		CCone<R> ZConeDual = ZCone.getDualizedCone_Parma().getSimplifiedConeFast();
		CMatrix<R> ZConeDualMatrix = DescriptorSetToMatrix(ZConeDual.getDescriptors(C));

		//ZConeDualMatrix.print();

		CMatrix<R> generatedM = ZConeDualMatrix * c;
		//generatedM.print();

		return CCone<R>(C, matrixToDescriptorSet(generatedM)).getSimplifiedConeFast();
	}
#endif

	template <class R> CCone<R> CCone<R>::getProjectionJones_ThroughPIP(std::vector<int> removedDims) {
		std::vector<int> retainedDims = getDimSetComplement(removedDims, this->getDimension());
		//CMatrix<R> (*U)(int dim) = &(CMatrix<R>::getUnitMatrix);
		//CMatrix<R> (*Z)(int rows, int cols) = &(CMatrix<R>::getZeroMatrix);
		CVector<R> (*UV)(int dim, int unitDim) = &(CVector<R>::getUnitVector);
		//CVector<R> (*CV)(int dim, R constant) = &(CVector<R>::getConstantVector);
		CVector<R> (*ZV)(int dim) = &(CVector<R>::getZeroVector);

		CMatrix<R> Y = matrix();
		CMatrix<R> d = CMatrix<R>(removedDims.size()); // Fill a matrix with the constraint coefficients corresponding to removed dimensions
		ITT(std::vector<int>, i, removedDims) d = d << Y.getColumn(*i);

		CMatrix<R> c = CMatrix<R>(retainedDims.size()); // Fill a matrix with the constraint coefficients corresponding to retained dimensions
		ITT(std::vector<int>, i, retainedDims) c = c << Y.getColumn(*i);


		//c.print();
		//d.print();

		// First get appropriate a vector
		// Unknown order: (g, a, u) // g = -t --> minimize
		/*CLpProblem<R> lpp = CLpProblem<R>(this->getDimension() + 1);
		CMatrix<R> lppCM = CMatrix<R>(this->getDimension() + 2);
		for (int q = 0; q < Y.getRowCount(); q++) {
			R sum = R::getZero();
			for (int r = 0; r < Y.getColumnCount(); r++) sum += Y[q][r] * Y[q][r];
			lppCM = lppCM >> (ZV(1) << CV(1, sum) << c[q] << -d[q]);
		}
		lppCM = lppCM >> (-UV(lppCM.getColumnCount(), 1)) // g <= 0
		              >> (UV(lppCM.getColumnCount(), 0) + UV(lppCM.getColumnCount(), 1)); // g >= -1

		lppCM.print();
		lpp.domConstrPolyheder = CFlat<CCone<R> >(1, CCone<R>(C, matrixToDescriptorSet(lppCM)));
		//lpp.maximize = true;
		CVector<typename R::RationalType>* lppSolution = lpp.solve();
		lppSolution->print();
		CVector<R> appropriateA = CVector<R>();
		R appropriateADenominatorMult = R::getOne();
		for (int q = 0; q < c.getColumnCount(); q++) {
			appropriateA.appendElement((*lppSolution)[1 + q].numerator);
			appropriateADenominatorMult = appropriateADenominatorMult * (*lppSolution)[1 + q].denominator;
		}
		delete lppSolution;
		appropriateA.print();*/

		// Now solve the actual PIP that will give us the resulting constraints
		/*c = (c << Z(c.getRowCount(), 1))
		    >> (appropriateA << -CV(1, appropriateADenominatorMult))
		    >> -(appropriateA << -CV(1, appropriateADenominatorMult));
		c.print();
		d = d >> ZV(d.getColumnCount()) >> ZV(d.getColumnCount());
		CMatrix<R> CT = c.getTransposedMatrix();
		CMatrix<R> DT = d.getTransposedMatrix();
		CT.print();
		DT.print();


		int dC = CT.getRowCount();
		int colsCT = CT.getColumnCount();
		int dD = DT.getRowCount();

		CMatrix<R> CM =   (CT << -U(dC))
		              >> -(CT << -U(dC))
		              >>  (DT << Z(dD, dC))
		              >> -(DT << Z(dD, dC))
		              >>  (U(colsCT) << Z(colsCT, dC));

		CM.print();

		CMatrix<R> CMO = ((CMatrix<R>) (UV(1, 0) << ZV(colsCT) << CVector<R>::getConstantVector(dC, -R::getOne())))
		                 >> -((CMatrix<R>) (UV(1, 0) << ZV(colsCT) << CVector<R>::getConstantVector(dC, -R::getOne())))
		                 >> (Z(CM.getRowCount(), 1) << CM);
		CMO = Z(CMO.getRowCount(), 1) << CMO
			//>> (CVector<R>::getConstantVector(1, -R::getOne()) << ZV(1 + consCount) << UV(dC, 0)) // New homog >= 1
			//>> -(CVector<R>::getConstantVector(1, -R::getOne()) << ZV(1 + consCount) << UV(dC, 0)) // New homog <= 1
			;

		CMO.print();

		CPipProblem<R> pp = CPipProblem<R>(colsCT + 1, dC);
		pp.integerSolution = false;
		pp.simplify = true;
		pp.domConstrFlat = CFlat<CModule<R> >(1, CModule<R>::universe(CMO.getColumnCount()));
		pp.domConstrPolyheder = CFlat<CCone<R> >(1, CCone<R>(C, matrixToDescriptorSet(CMO)));

		CPipQuast<R>* solution = pp.solve();

		ASSERT(solution != NULL);
		if (solution != NULL) solution->print();
		*/

		// calc [C b]a == c a
		/*CMatrix<R> CM2 = CMatrix<R>(2 + d.getColumnCount() + c.getColumnCount());
		for (int q = 0; q < c.getRowCount(); q++) {
			CVector<R> v = CVector<R>();
			R sum = R::getZero();
			for (int r = 0; r < appropriateA.getLength(); r++) {
				sum = sum + (c[q][r] * appropriateA[r]);
			}
			v.appendElement(R::getZero());
			v.appendElement(-sum);
			v = v << -(appropriateADenominatorMult*d[q]) << (appropriateADenominatorMult*c[q]);
			CM2 = CM2 >> v;
		}

		CM2.print();
		CPipProblem<R> pp = CPipProblem<R>(1 + d.getColumnCount(), c.getColumnCount());
		pp.integerSolution = false;
		pp.simplify = true;
		//pp.domConstrFlat = CFlat<CModule<R> >(1, CModule<R>::universe(CM2.getColumnCount()));
		pp.domConstrPolyheder = CFlat<CCone<R> >(1, CCone<R>(C, matrixToDescriptorSet(CM2)));

		CPipQuast<R>* solution = pp.solve();

		ASSERT(solution != NULL);
		if (solution != NULL) solution->print();*/

		// calc [C b]a == c a
		CMatrix<R> CM2 = CMatrix<R>(2 + d.getColumnCount() + c.getColumnCount());
		for (int q = 0; q < c.getRowCount(); q++) {
			/*CVector<R> v = CVector<R>();
			R sum = R::getZero();
			for (int r = 0; r < appropriateA.getLength(); r++) {
				sum = sum + (c[q][r] * appropriateA[r]);
			}*/
			CM2 = CM2 >> (UV(2, 1) << (d[q]) << (c[q]));
			CM2 = CM2 >> (ZV(2) << (d[q]) << (c[q]));
		}
		CM2 = CM2 >> -(UV(CM2.getColumnCount(), 1));

		//CM2.print();
		CPipProblem<R> pp = CPipProblem<R>(1 + d.getColumnCount(), c.getColumnCount());
		pp.integerSolution = true;
		pp.simplify = true;
		//pp.domConstrFlat = CFlat<CModule<R> >(1, CModule<R>::universe(CM2.getColumnCount()));
		pp.domConstrPolyheder = CFlat<CCone<R> >(1, CCone<R>(C, matrixToDescriptorSet(CM2)));

		CPipQuast<R>* solution = pp.solve();

		ASSERT(solution != NULL);
		//if (solution != NULL) solution->print();

		CMatrix<R> newCs = extractConstraintMatrixFromPipQuast(retainedDims.size(), solution);
		CCone<R> result = CCone<R>(C, matrixToDescriptorSet(newCs)).getSimplifiedConeFast();
		//result.print();
		return result;
	}

	template <class R> CCone<R> CCone<R>::getProjectionFourierMotzkin(int projectionDim, bool preSimplify) {
		ASSERT((projectionDim >= 0) && (projectionDim < this->getDimension()));

		CMatrix<R> Y = preSimplify ? getSimplifiedConeFast().matrix() : matrix();
		const CMatrix<R>& X = Y;
		//X.print();

		CMatrix<R> result = CMatrix<R>(this->getDimension() - 1);
		vector<int> posList, negList; // Lists with constraints with zero-coeff, pos coeff, neg coeff for projectionDim respectively
		for (int r = 0; r < X.getRowCount(); r++) { R pivot = X(r, projectionDim);
			if (pivot == 0) {
				result.addRow(X[r].getSubVector(0, projectionDim)
					                 << X[r].getSubVector(projectionDim + 1, result.getColumnCount() - projectionDim));
			} else if (pivot > 0) { posList.push_back(r); }
			else { negList.push_back(r); }
		}

		//result.print();
		int firstNewConstraint = result.getRowCount();

		// Recombine pos & neg
		//CInteger(posList.size()*negList.size());
		for (unsigned int p = 0; p < posList.size(); p++) { int rP = posList[p];
			for (unsigned int n = 0; n < negList.size(); n++) { int rN = negList[n];
				CVector<R>& row = *(new CVector<R>(result.getColumnCount()));
				bool isZero = true;
				for (int c = 0; c < result.getColumnCount(); c++) { int d = (c < projectionDim) ? c : c + 1;
					row[c] = (X(rP, projectionDim) * X(rN, d) - X(rP, d) * X(rN, projectionDim));
					isZero = isZero && (row[c] == 0);
				}
				if (!isZero) {
					result.addRow(&row);
				} else delete &row;
			}
		}

		//result.print();

		return CCone<R>(C, matrixToDescriptorSet(result)).getSimplifiedConeFast_Ex(firstNewConstraint);
	}

#ifdef PPL_LINK
	template <class G> CCone<G> CCone<G>::getDualizedCone_Parma() {
		Parma_Polyhedra_Library::Constraint_System pcm = (Parma_Polyhedra_Library::Constraint_System) matrix();

		Parma_Polyhedra_Library::C_Polyhedron poly = Parma_Polyhedra_Library::C_Polyhedron(pcm);


		//poly.update_generators();

		Parma_Polyhedra_Library::Generator_System gensys = poly.generators();

		return CCone<G>(C, matrixToDescriptorSet(CMatrix<G>::fromParmaGeneratorSystem(gensys)));
	}
#endif

	template <class G> CCone<G> CCone<G>::getDualizedCone() {
		CMatrix<G> (*U)(int dim) = &(CMatrix<G>::getUnitMatrix);
		CMatrix<G> (*Z)(int rows, int cols) = &(CMatrix<G>::getZeroMatrix);

		CMatrix<G> consT = getSimplifiedConeFast().matrix().getTransposedMatrix();
		//constraintMatrix.print();
		CCone<G> w = CCone<G>(C, matrixToDescriptorSet(
				(-consT <<  U(this->getDimension()))
			 >> ( consT << -U(this->getDimension()))
	         >> ( U(consT.getColumnCount())      <<  Z(consT.getColumnCount(), this->getDimension()))));

		//w.print();
		printf("[%i|", w.matrix().getColumnCount());
		for (int q = 0; q < consT.getColumnCount(); q++) {
			w = w.getProjectionFourierMotzkin(0, (q == 0));
			printf(":%i", w.matrix().getRowCount());
			//w.print();
		}
		printf(":|%i]", w.matrix().getColumnCount());

		//w.print();
		//w = w.getProjectionJones_ThroughPIP(CCone<G>::getInterval(0, consT.getColumnCount()));
		//w.print();

		return w;
	}


	/* Simplification procedure:
	 * (Zero constraints are redundant)
	 * - Constraints which only differ in (positive) scale can be represented by one of them
	 * - LP for each constraint to check redundance in retainNonRedundant()
	 */

	template <class R> bool pVectorCompare(CVector<R>* v1, CVector<R>* v2) {
		return *v2 < *v1;
	}

	template <class R> bool pVectorRedundanceStructCompare(const VectorRedundanceStruct<R>& vrs1, const VectorRedundanceStruct<R>& vrs2) {
		return (*(vrs1.v) < *(vrs2.v))
		       ? true
		       : ((*(vrs1.v) == *(vrs2.v)) ? (vrs1.certainlyNonRedundant && !vrs2.certainlyNonRedundant) : false);
	}

	template <class R> CCone<R> CCone<R>::getSimplifiedConeFast_Ex(int numberOfCertainlyNonRedundant) {
		CMatrix<R> consMatrix = matrix();

		vector<VectorRedundanceStruct<R> > constrs = vector<VectorRedundanceStruct<R> >();

		// First divide each constraint by the gcd
		for (int q = 0; q < consMatrix.getRowCount(); q++) {
			consMatrix[q].divideByGCDOfElements();
			constrs.push_back(VectorRedundanceStruct<R>(&consMatrix[q], q < numberOfCertainlyNonRedundant));
		}

		std::sort(constrs.begin(), constrs.end(), pVectorRedundanceStructCompare<R>);

		// Build a list of less redundant constraints
		CVector<R>* lastV = NULL;
		list<VectorRedundanceStruct<R>* > remConstrs = list<VectorRedundanceStruct<R>* >();
		for (unsigned int q = 0; q < constrs.size(); q++) {
			if ((q == 0) || (*(constrs[q].v) != *lastV)) {
				remConstrs.push_back(&constrs[q]);
				lastV = constrs[q].v;
			}
		}

		// Remove redundant constraints using integer programming
		remConstrs = retainNonRedundant(remConstrs);

		CMatrix<R> resultConstr = CMatrix<R>(this->getDimension());
		ITT(list<VectorRedundanceStruct<R>* >, lit, remConstrs) resultConstr.addRow(*((*lit)->v));

		return CCone<R>(C, matrixToDescriptorSet(resultConstr));
	}

	template <class R> CCone<R> CCone<R>::getSimplifiedConeFast() {
		CMatrix<R> consMatrix = matrix();

		vector<CVector<R>* > constrs = vector<CVector<R>* >();

		// First divide each constraint by the gcd
		for (int q = 0; q < consMatrix.getRowCount(); q++) {
			consMatrix[q].divideByGCDOfElements();
			constrs.push_back(&consMatrix[q]);
		}

		std::sort(constrs.begin(), constrs.end(), pVectorCompare<R>);

		// Build a list of less redundant constraints
		CVector<R>* lastV = NULL;
		list<CVector<R>* > remConstrs = list<CVector<R>* >();
		for (unsigned int q = 0; q < constrs.size(); q++) {
			if ((q == 0) || (*constrs[q] != *lastV)) {
				remConstrs.push_back(constrs[q]);
				lastV = constrs[q];
			}
		}

		// Remove redundant constraints using integer programming
		remConstrs = retainNonRedundant(remConstrs);

		CMatrix<R> resultConstr = CMatrix<R>(this->getDimension());
		typedef typename list<CVector<R>* >::iterator listIt;
		for (listIt lit = remConstrs.begin(); lit != remConstrs.end(); ++lit) resultConstr.addRow(**lit);

		return CCone<R>(C, matrixToDescriptorSet(resultConstr));
	}

	template <class R>
	list<CVector<R>* > CCone<R>::retainNonRedundant(list<CVector<R>* > cons) const {
		typedef typename list<CVector<R>* >::iterator listIt;

		list<CVector<R>* > nonRed = list<CVector<R>* >();

		while (cons.size() > 0) {
			list<CVector<R>* > othC;

			listIt np = cons.begin();
			CVector<R>* currC = *np;
			for (listIt slit = nonRed.begin(); slit != nonRed.end(); ++slit) othC.push_back(*slit);
			for (listIt slit = ++np; slit != cons.end(); ++slit) othC.push_back(*slit);
			if (!isRedundant(*currC, othC)) nonRed.push_back(currC);
			cons.erase(cons.begin());
		}

		return nonRed;
	}

	template <class R>
	list<VectorRedundanceStruct<R>* > CCone<R>::retainNonRedundant(list<VectorRedundanceStruct<R>* > cons) const {
		typedef typename list<VectorRedundanceStruct<R>* >::iterator listIt;

		list<VectorRedundanceStruct<R>* > nonRed = list<VectorRedundanceStruct<R>* >();

		while (cons.size() > 0) {
			list<CVector<R>* > othC;

			listIt np = cons.begin();
			VectorRedundanceStruct<R>* currC = *np;
			if (currC->certainlyNonRedundant) {
				nonRed.push_back(currC);
			} else {
				for (listIt slit = nonRed.begin(); slit != nonRed.end(); ++slit) othC.push_back((*slit)->v);
				for (listIt slit = ++np; slit != cons.end(); ++slit) othC.push_back((*slit)->v);
				if (!isRedundant(*(currC->v), othC)) nonRed.push_back(currC);
			}
			cons.erase(cons.begin());
		}

		return nonRed;
	}

	template <class R>
	bool CCone<R>::isRedundant(CVector<R>& c, list<CVector<R>* > otherC) const {
		typedef typename R::RationalType Rational;
		typedef typename list<CVector<R>* >::iterator listIt;

		CVector<R> (*UV)(int dim, int unitDim) = &(CVector<R>::getUnitVector);
		CVector<R> (*ZV)(int dim) = &(CVector<R>::getZeroVector);

		if (otherC.size() == 0) return (c == ZV(c.getLength()));

		CLpProblem<R> ipp = CLpProblem<R>(c.getLength() + 1);

		CMatrix<R> domConM = CMatrix<R>(c.getLength() + 2);
		for (listIt lit = otherC.begin(); lit != otherC.end(); ++lit) domConM.addRow(ZV(2) << (**lit));
		domConM.addRow(UV(2, 1) << -c); // h >= c.x~ ( >= 0 )
		domConM.addRow(UV(1, 0) << UV(1, 0) << ZV(c.getLength())); // h >= -1   XXX:


		//domConM.print(); //c.print();

		ipp.setDomainPolyheder(CFlat<CCone<R> >(1, CCone<R>(C, matrixToDescriptorSet(domConM))));

		CVector<Rational>* solution = ipp.solve();

		ASSERT(solution != NULL);

		/*ParmaMIPProblem plpp = ParmaMIPProblem(domConM.getColumnCount() - 1, (ParmaConstraintSystem) domConM,
				                (ParmaLinearExpression) UV(c.getLength() + 1, 1),
				                Parma_Polyhedra_Library::MINIMIZATION);

		plpp.solve();*/

		//print();

		//bool redundant = CVector<R>::fromParmaGenerator(plpp.optimizing_point())[0].isZero();
		//(*solution).print();
		bool redundant = ((*solution)[0] == Rational::getZero());
		delete solution;


		return redundant;
	}

	template <class R> CCone<R> CCone<R>::getReflexiveTransitiveClosure(int spaceDim) {
		int affineness = this->getDimension() - 2*spaceDim;
		CMatrix<R> (*Z)(int rows, int cols) = &(CMatrix<R>::getZeroMatrix);
		CMatrix<R> (*U)(int dim) = &(CMatrix<R>::getUnitMatrix);
		//CVector<R> (*UV)(int dim, int unitDim) = &(CVector<R>::getUnitVector);
		//CVector<R> (*ZV)(int dim) = &(CVector<R>::getZeroVector);

		CMatrix<R> consM = matrix();
		CCone<R> c = CCone<R>::fromConstraints((consM << Z(consM.getRowCount(), spaceDim))
				  >>  (Z(spaceDim, affineness) << -U(spaceDim) << U(spaceDim) << -U(spaceDim))
				  >> -(Z(spaceDim, affineness) << -U(spaceDim) << U(spaceDim) << -U(spaceDim)));

		for (int q = 0; q < 2*spaceDim; q++) c = c.getProjectionFourierMotzkin(affineness, q == 0);
		for (int q = 0; q < affineness; q++) c = c.getProjectionFourierMotzkin(0, false);

		return c;
	}

	template <class R>
	const CMatrix<R> CCone<R>::matrix() {
		return DescriptorSetToMatrix(CHadron<CUnidirectionalVectorDescriptorSet<R> >::getDescriptors(C));
	}

	/*template <class R>
	std::vector<int> CCone<R>::getDimSetComplement(std::vector<int> dims, int dimension) {
		std::vector<bool> complDim = std::vector<bool>(dimension, true);
		for (unsigned int q = 0; q < dims.size(); q++) complDim[dims[q]] = false;
		std::vector<int> resultDims;
		for (unsigned int q = 0; q < complDim.size(); q++) if (complDim[q]) resultDims.push_back(q);

		return resultDims;
	}*/

	/*
	template <class R>
	std::vector<int> CCone<R>::getInterval(int start, int length) {
		std::vector<int> resultDims;

		for (int q = 0; q < length; q++) resultDims.push_back(start + q);

		return resultDims;
	}
	*/


}

#endif /*CONE_H_*/
