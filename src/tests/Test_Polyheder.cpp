
#include "../cute/cute.h"
#include "../cute/ide_listener.h"

#include "Test_Funcs.h"

#include "../partyql/basicmath/scalar/Rational.h"
#include "../partyql/basicmath/scalar/Interval.h"
//#include "../partyql/basicmath/Polyheder.h"

#include "../partyql/core/util/Time.h"


using namespace AlgoTrans;


namespace Test_LimsLinearisation {
	typedef CRational<I> IRat;
	typedef CInterval<IRat> IRatIntv;

	void A() {
		Polyheder pIn = IP(C, 1, IV_(2, -7,  1));
		Polyheder exp = Polyheder::universe(1); // The entire space

		Polyheder affineHull = pIn.getBidirectionalHull();
		Polyheder affineHull_Lim = pIn.getBidirectionalHull_Lim();

		ASSERT(affineHull == exp);
		ASSERT(affineHull_Lim == exp);
	}

	void B() {
		Polyheder pIn = IP(C, 2, IV_(2, -7,  1),
							            IV_(2,  7,  1));
		Polyheder exp = Polyheder::universe(1); // The entire space

		Polyheder affineHull = pIn.getBidirectionalHull();
		Polyheder affineHull_Lim = pIn.getBidirectionalHull_Lim();

		pIn.print();
		affineHull.print();
		ASSERT(affineHull == exp);
		ASSERT_EQ(affineHull_Lim, exp);
	}

	void B2() {
		Polyheder pIn = IP(C, 2, IV_(2, -7,  1),
							            IV_(2,  7, -1));
		Polyheder exp = Polyheder::universe(1);

		Polyheder affineHull = pIn.getBidirectionalHull();
		Polyheder affineHull_Lim = pIn.getBidirectionalHull_Lim();

		ASSERT(affineHull == exp);
		ASSERT(affineHull_Lim == exp);
	}

	void C_() {
		Polyheder pIn = IP(C, false, 4, IV_(3,  2,  1, -1),  // x - y >= -2
			                            IV_(3,  2, -1,  1),  // x - y <=  2
				                        IV_(3,  8, -1, -1),  // x + y <=  8
				                        IV_(3, -4,  1,  1)); // x + y >=  4

		Polyheder exp = Polyheder::universe(2); // The entire space

		Polyheder affineHull = pIn.getBidirectionalHull();
		Polyheder affineHull_Lim = pIn.getBidirectionalHull_Lim();

		ASSERT(affineHull == exp);
		ASSERT(affineHull_Lim == exp);
	}

	void D() {
		Polyheder pIn = IP(C, false, 5, IV_(3,  2,  1, -1),  // x - y >= -2
			                            IV_(3,  2, -1,  1),  // x - y <=  2
				                        IV_(3,  8, -1, -1),  // x + y <=  8
				                        IV_(3, -4,  1,  1),  // x + y >=  4
				                        IV_(3, -8,  1,  1)); // x + y >=  8

		Polyheder exp = IP(C, true, 1, IV_(3, -8,  1,  1));

		Polyheder affineHull = pIn.getBidirectionalHull();
		Polyheder affineHull_Lim = pIn.getBidirectionalHull_Lim();

		ASSERT(affineHull == exp);
		ASSERT(affineHull_Lim == exp);
	}

	/* void E() { // Empty poly -- does not work with Lim's Linearisation due to the use of Farkas' Lemma
		IPolyheder pIn = IP(3, IV_(3,   9,  1,  0),  // x >=  -9
			                   IV_(3,  -9, -1,  0),  // x <=  -9
				               IV_(3, -10, -1,  0)); // x <= -10

		IModule affineHull = pIn.getAffineHull();

		IModule exp = IMod_Con(1, IV_(3, -1, 0, 0));

		ASSERT(affineHull == exp);
	} */

	void F() {
		Polyheder pIn = IP(C, false, 2, IV_(3,   9,  1,  0),  // x >=  -9
			                            IV_(3,  -9, -1,  0)); // x <=  -9

		Polyheder exp = IP(C, true, 1, IV_(3,  9,  1,  0));

		Polyheder affineHull = pIn.getBidirectionalHull();
		Polyheder affineHull_Lim = pIn.getBidirectionalHull_Lim();

		ASSERT(affineHull == exp);
		ASSERT(affineHull_Lim == exp);
	}

	void hyperCubes() {
		IMatrix (*U)(int dim) = &(CMatrix<I>::getUnitMatrix);
		IMatrix (*Z)(int rows, int cols) = &(CMatrix<I>::getZeroMatrix);
		IVector (*UV)(int dim, int unitDim) = &(CVector<I>::getUnitVector);

		int maxDim = 3;
		for (int q = 1; q < maxDim; q++) {
			printf("%i: ", q);
			CMatrix<I> cons = Z(q, 1) << U(q); // x >= 0
			for (int r = 0; r < q; r++) cons.addRow(UV(1, 0) << (-UV(q, r)));

			IPolyheder p = IPolyheder(cons);

			double tPrev, tDuration;
			IFlatModule *affineHull, *affineHullLim;
			TIME( affineHullLim = new IFlatModule(p.getAffineHull_LimsLinearisation()); , "AffineHull_Lim" );
			TIME( affineHull = new IFlatModule(p.getAffineHull()); , "AffineHull" );

			IFlatModule exp = IFMod_Con(1, 0, q + 1);

			//ASSERT(affineHullLim == exp);
			ASSERT(*affineHull == exp);

			delete affineHull;
			delete affineHullLim;
		}
	}

	void hyperPyramids() {
		IMatrix (*U)(int dim) = &(CMatrix<I>::getUnitMatrix);
		IMatrix (*Z)(int rows, int cols) = &(CMatrix<I>::getZeroMatrix);
		IVector (*UV)(int dim, int unitDim) = &(CVector<I>::getUnitVector);
		IVector (*CV)(int dim, I constant) = &(CVector<I>::getConstantVector);

		int maxDim = 3;
		for (int q = 1; q < maxDim; q++) {
			printf("%i: ", q);
			CMatrix<I> cons = Z(q, 1) << U(q); // x >= 0
			cons.addRow(UV(1, 0) << -CV(q, 1));

			IPolyheder p = IPolyheder(cons);

			double tPrev, tDuration;
			IFlatModule *affineHull, *affineHullLim;
			TIME( affineHullLim = new IFlatModule(p.getAffineHull_LimsLinearisation()); , "AffineHull_Lim" );
			TIME( affineHull = new IFlatModule(p.getAffineHull()); , "AffineHull" );

			IFlatModule exp = IFMod_Con(1, 0, q + 1);

			//ASSERT(affineHullLim == exp);
			ASSERT(*affineHull == exp);

			delete affineHull;
			delete affineHullLim;
		}
	}

	/*void G() { // Empty poly -- does not work with Lim's Linearisation due to the use of Farkas' Lemma
		IPolyheder pIn = IP(2, IV_(3,   9,  1,  0),  // x >=  -9
			                   IV_(3, -10, -1,  0)); // x <= -10

		IModule affineHull = pIn.getAffineHull();

		ASSERT(affineHull.isEmpty());
	}

	void H() { // Empty poly -- does not work with Lim's Linearisation due to the use of Farkas' Lemma
		IPolyheder pIn = IP(2, IV_(2,   9,  1),  // x >=  -9
			                   IV_(2, -10, -1)); // x <= -10

		(!pIn.getAffineHull()).print();

		IModule affineHull = pIn.getAffineHull();

		ASSERT(affineHull.isEmpty());
	}*/
}

namespace Test_Polyheder {
	void equality_A() {
		ASSERT(IP(C, false, 1, IV_(2, 5, 1)) == IP(C, false, 1, IV_(2, 5, 1)));

		ASSERT(IP(C, false, 2, IV_(3,  9,  1,  0),
		             IV_(3, -9, -1,  0)) == IP(C, false, 2, IV_(3,  9,  1,  0),
			                                      IV_(3, -9, -1,  0)));

		ASSERT(IP(C, false, 4, IV_(3,  9,  1,  0),
  		             IV_(3, -9, -1,  0),
		             IV_(3,  9,  0,  1),
		             IV_(3, -9,  0, -1)) == IP(C, false, 4, IV_(3,  9,  1,  0),
		            		 					  IV_(3, -9, -1,  0),
		            		 					  IV_(3,  9,  0,  1),
		            		 					  IV_(3, -9,  0, -1)));
	}

	void inequality_A() {
		ASSERT(IP(C, false, 1, IV_(2, 5, 1)) != IP(C, false, 1, IV_(2, 5, 2)));
		ASSERT(IP(C, false, 1, IV_(2, 5, 1)) != IP(C, false, 1, IV_(2, 6, 1)));
		ASSERT(IP(C, false, 1, IV_(2, 5, 1)) != IP(C, false, 1, IV_(2, -5, 1)));
		ASSERT(IP(C, false, 1, IV_(2, 5, 1)) != IP(C, false, 1, IV_(2, 5, -1)));

		ASSERT(IP(C, false, 2, IV_(3,  9,  1,  0),
		             IV_(3, -9, -1,  0)) != IP(C, false, 2, IV_(3, 10,  1,  0),
			                                      IV_(3, -9, -1,  0)));

		ASSERT(IP(C, false, 2, IV_(3, 10,  1,  0),
		             IV_(3, -9, -1,  0)) != IP(C, false, 2, IV_(3,  9,  1,  0),
			                                      IV_(3, -9, -1,  0)));

		ASSERT(IP(C, false, 2, IV_(3,  9,  1,  0),
		             IV_(3, -9, -1,  0)) != IP(C, false, 2, IV_(3,  9,  2,  0),
			                                      IV_(3, -9, -1,  0)));

		ASSERT(IP(C, false, 2, IV_(3,  9,  1,  0),
		             IV_(3, -9, -1,  0)) != IP(C, false, 2, IV_(3,  9, -1,  0),
			                                      IV_(3, -9, -1,  0)));

		ASSERT(IP(C, false, 2, IV_(3,  9,  1,  0),
		             IV_(3, -9, -1,  0)) != IP(C, false, 2, IV_(3,  9,  1,  1),
			                                      IV_(3, -9, -1,  0)));

		ASSERT(IP(C, false, 4, IV_(3,  9,  1,  0),
  		             IV_(3, -9, -1,  0),
		             IV_(3,  9,  0,  1),
		             IV_(3, -9,  0, -1)) != IP(C, false, 3, IV_(3,  9,  1,  0),
		            		 					  IV_(3, -9, -1,  0),
		            		 					  IV_(3,  9,  0,  1)));
	}
}

namespace Test_ConeDual {
	void A() {
		Polyheder a = IP(C, false, 1, IV_(2, 5, -1));
		Polyheder a_d = a.getDual();

		Polyheder a_exp = IP(C, false, 2, IV_(2, 1, 5),
				                 IV_(2, 0, -1));

		ASSERT(a_d == a_exp);


		Polyheder b = IP(C, false, 1, IV_(2, 5, 1));
		Polyheder b_d = b.getDual();

		Polyheder b_exp = IP(C, false, 2, IV_(2, 1, -5),
				                 IV_(2, 0, 1));

		ASSERT(b_d == b_exp);


		Polyheder c = IP(C, false, 1, IV_(2, -5, -1));
		Polyheder c_d = c.getDual();

		Polyheder c_exp = IP(C, false, 2, IV_(2, 1, -5),
				                 IV_(2, 0, -1));

		ASSERT(c_d == c_exp);


		Polyheder d = IP(C, false, 1, IV_(2, -5, 1));
		Polyheder d_d = d.getDual();

		Polyheder d_exp = IP(C, false, 2, IV_(2, 1, 5),
				                 IV_(2, 0, 1));

		ASSERT(d_d == d_exp);
	}

	void polyheder_intersection() {
		Polyheder a = IP(C, false, 1, IV_(2, 5, 1)); // x >= -5
		Polyheder b = IP(C, false, 1, IV_(2, 5, -1)); // x <= 5
		Polyheder abIntsect = IP(C, false, 2, IV_(2, 5, 1), IV_(2, 5, -1)); // x >= -5 && x <= 5

		ASSERT((a && b) == abIntsect);
	}

	void polyheder_convexHull() {
		Polyheder a = IP(C, false, 2, IV_(2, 5, 1),    // x >= -5
				             IV_(2, -3, -1)); // x <= -3
		Polyheder b = IP(C, false, 2, IV_(2, 5, -1),  // x <= 5
				             IV_(2, -3, 1)); // x >= 3
		Polyheder abConvexHull = IP(C, false, 2, IV_(2, 5, 1), IV_(2, 5, -1)); // x >= -5 && x <= 5

		ASSERT((a + b) == abConvexHull);
	}
}

namespace Test_FlatConcatenation {
	void flatConcatenation() {
		//: DS<D<I>>(3, D<I>( 1 1 0 == ), D<I>( 0 0 1 == ))
		//: DS<D<I>>(2, D<I>( -1  0 == ), D<I>(  0 -1 == ))
	}
}

namespace Test_ConvexHull {
	void hyperCubes() {
		IMatrix (*U)(int dim) = &(CMatrix<I>::getUnitMatrix);
		IMatrix (*Z)(int rows, int cols) = &(CMatrix<I>::getZeroMatrix);
		IVector (*UV)(int dim, int unitDim) = &(CVector<I>::getUnitVector);

		int maxDim = 3;
		for (int q = 1; q < maxDim; q++) {
			printf("%i: ", q);
			CMatrix<I> cons = Z(q, 1) << U(q); // x >= 0
			for (int r = 0; r < q; r++) cons.addRow(UV(1, 0) << (-UV(q, r)));

			IPolyheder p = IPolyheder(cons);

			double tPrev, tDuration;
			IPolyheder *convexHull;
			TIME( convexHull = new IPolyheder(p + p); , "ConvexHull" );

			IPolyheder exp = p;

			//ASSERT(affineHullLim == exp);
			ASSERT(*convexHull == exp);

			delete convexHull;
		}
	}

	/*void hyperPyramids() {
		IMatrix (*U)(int dim) = &(CMatrix<I>::getUnitMatrix);
		IMatrix (*Z)(int rows, int cols) = &(CMatrix<I>::getZeroMatrix);
		IVector (*UV)(int dim, int unitDim) = &(CVector<I>::getUnitVector);
		IVector (*CV)(int dim, I constant) = &(CVector<I>::getConstantVector);

		int maxDim = 7;
		for (int q = 1; q < maxDim; q++) {
			printf("%i: ", q);
			CMatrix<I> cons = Z(q, 1) << U(q); // x >= 0
			cons.addRow(UV(1, 0) << -CV(q, 1));

			IPolyheder p = IPolyheder(cons);

			double tPrev, tDuration;
			IFlatModule *affineHull, *affineHullLim;
			TIME( affineHullLim = new IFlatModule(p.getAffineHull_LimsLinearisation()); , "AffineHull_Lim" );
			TIME( affineHull = new IFlatModule(p.getAffineHull()); , "AffineHull" );

			IFlatModule exp = IFMod_Con(1, 0, q + 1);

			//ASSERT(affineHullLim == exp);
			ASSERT(*affineHull == exp);

			delete affineHull;
			delete affineHullLim;
		}
	}*/
}

cute::suite* Test_LimsLinearisation_runSuite(){
	AlgoTrans::timeConvexHull_J = 0.0;
	AlgoTrans::timeConvexHull_FM = 0;

	cute::suite& s = *(new cute::suite("Polyheder"));

	s.push_back(CUTE(Test_LimsLinearisation::A));
	s.push_back(CUTE(Test_LimsLinearisation::B));
	s.push_back(CUTE(Test_LimsLinearisation::B2));
	s.push_back(CUTE(Test_LimsLinearisation::C_));
	s.push_back(CUTE(Test_LimsLinearisation::D));
	s.push_back(CUTE(Test_LimsLinearisation::F));
	s.push_back(CUTE(Test_LimsLinearisation::hyperCubes));
	s.push_back(CUTE(Test_LimsLinearisation::hyperPyramids));

	s.push_back(CUTE(Test_Polyheder::equality_A));
	s.push_back(CUTE(Test_Polyheder::inequality_A));

	s.push_back(CUTE(Test_ConeDual::A));
	s.push_back(CUTE(Test_ConeDual::polyheder_intersection));
	s.push_back(CUTE(Test_ConeDual::polyheder_convexHull));

//	s.push_back(CUTE(Test_ConvexHull::hyperCubes));

	s.push_back(CUTE(Test_FlatConcatenation::flatConcatenation));
	//s.push_back(CUTE(Test_ConvexHull::hyperPyramids));

	// Skipped -- Empty polyhedra not processable by Lim's Linearisation
	//s.push_back(CUTE(Test_LimsLinearisation::H));
	//s.push_back(CUTE(Test_LimsLinearisation::G));
	//s.push_back(CUTE(Test_LimsLinearisation::E));

	return &s;
}


