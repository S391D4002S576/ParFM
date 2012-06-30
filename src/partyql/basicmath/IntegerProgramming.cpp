#ifdef INTEGERPROGRAMMING2_H_

namespace AlgoTrans {
	template <class R> PipMatrix*
	CPipProblem<R>::getPipConstrMatrixFromAffineSetAndPolyheder(int colCount, const CMatrix<R>& eqs, const CMatrix<R>& ineqs) {
	    //ASSERT(eqs.getColumnCount() == ineqs.getColumnCount());
		//eqs.print();
		//ineqs.print();

		PipMatrix* result = pip_matrix_alloc(eqs.getRowCount() + ineqs.getRowCount(), colCount + 1);
	    Entier* row;

	    for (int r = 0; r < eqs.getRowCount(); r++) {
	        row = result->p[r];
	        value_assign(row[0], CInteger(0).getGMPData());
	        for (int c = 0; c < eqs.getColumnCount(); c++) value_assign(row[1 + c], eqs[r][c].getGMPData());
	    }

	    for (int r = 0; r < ineqs.getRowCount(); r++) {
	        row = result->p[eqs.getRowCount() + r];
	        value_assign(row[0], CInteger(1).getGMPData());
	        for (int c = 0; c < ineqs.getColumnCount(); c++) value_assign(row[1 + c], ineqs[r][c].getGMPData());
	    }

	    //pip_matrix_print(stdout, result);

	    return result;
	}

	template <class R> CPipQuast<R>*
	CPipProblem<R>::solve() {
		//ASSERT(paramConstrFlat.getSpaceDimension() == paramConstrPolyheder.getSpaceDimension());
		//ASSERT(paramConstrFlat.getSpaceDimension() + 1 == domConstrFlat.getAffineness());
		//ASSERT(domConstrFlat.getSpaceDimension() == domConstrPolyheder.getSpaceDimension());

		int aff = parameterCount + 1;
		int sd = unknownCount;
		CMatrix<R> domCFM = DescriptorSetToMatrix(domConstrFlat.getDescriptors(C));
		CMatrix<R> pipDCFM = (domCFM.getSubMatrixColumns(1, sd + aff - 1) << domCFM.getSubMatrixColumns(0,1));
		CMatrix<R> domCPM = DescriptorSetToMatrix(domConstrPolyheder.getDescriptors(C));
		CMatrix<R> pipDCPM = (domCPM.getSubMatrixColumns(1, sd + aff - 1) << domCPM.getSubMatrixColumns(0,1));

		PipMatrix* domConstrMatrix = getPipConstrMatrixFromAffineSetAndPolyheder(sd + aff, pipDCFM, pipDCPM);

		CMatrix<R> parCFM = DescriptorSetToMatrix(paramConstrFlat.getDescriptors(C));
		CMatrix<R> pipPCFM = (parCFM.getSubMatrixColumns(1, aff - 1) << parCFM.getSubMatrixColumns(0, 1));
		CMatrix<R> parCPM = DescriptorSetToMatrix(paramConstrPolyheder.getDescriptors(C));
		CMatrix<R> pipPCPM = (parCPM.getSubMatrixColumns(1, aff - 1) << parCPM.getSubMatrixColumns(0, 1));

		PipMatrix* parConstrMatrix = getPipConstrMatrixFromAffineSetAndPolyheder(aff, pipPCFM, pipPCPM);

		// Fill options
		PipOptions* options = pip_options_init();
		options->Nq = integerSolution ? 1 : 0;
		options->Verbose = -1;
		options->Simplify = simplify ? 1 : 0;
		options->Deepest_cut = deepestCut ? 1 : 0;
		options->Maximize = maximize ? 1 : 0;
		options->Urs_parms = positiveParameters ? 0 : -1;
		options->Urs_unknowns = positiveUnknowns ? 0 : -1;

		// Solve
		PipQuast* solution = pip_solve(domConstrMatrix, parConstrMatrix, bigParameter, options);

		// Free matrices and options
		pip_options_free(options);
		pip_matrix_free(parConstrMatrix);
		pip_matrix_free(domConstrMatrix);

		// Process solution
		CPipQuast<R>* resultQuast = CPipQuast<R>::fromPipQuast(solution);

		pip_quast_free(solution);
		pip_close();

		return resultQuast;
	}

	template <class R>
	CPipProblem<R>::CPipProblem(int iUnknownCount, int iParameterCount)
	: unknownCount(iUnknownCount), parameterCount(iParameterCount),
	  domConstrFlat(CFlat<CModule<R> >::universe(iParameterCount + 1, iUnknownCount)),
	  domConstrPolyheder(CFlat<CCone<R> >::universe(iParameterCount + 1, iUnknownCount)),
	  paramConstrFlat(CFlat<CModule<R> >::universe(1, iParameterCount)),
	  paramConstrPolyheder(CFlat<CCone<R> >::universe(1, iParameterCount))
	{
		integerSolution = true;
		simplify = true;
		deepestCut = true;
		maximize = false;
		positiveParameters = false;
		positiveUnknowns = false;
		bigParameter = -1;
	}

	template <class R>
	CVector<typename R::RationalType>* CIpProblem<R>::solve() {
		typedef typename R::RationalType Rational;

		CPipQuast<R>* solution = CPipProblem<R>::solve();
		// XXX: Memory leak

		//solution->print();
		//ASSERT(solution->condition == NULL);

		CMatrix<Rational>* solM = solution->getBigParmSolution(0);
		//solM->print();
		if (solM == NULL) {
			delete solution;

			return NULL;
		} else {
			CMatrix<Rational> solC = solM->getSubMatrixColumns(solM->getColumnCount() - 1, 1);
			ASSERT(solC.getColumnCount() == 1);

			CVector<Rational>* result = new CVector<Rational>(solC.getRowCount());

			for (int q = 0; q < solM->getRowCount(); q++) (*result)[q] = solC[q][0];

			delete solution;

			return result;
		}
	}

	template <class R>
	CVector<typename R::RationalType>* CLpProblem<R>::solve() {
		typedef typename R::RationalType Rational;

		CPipQuast<R>* solution = CPipProblem<R>::solve();
		// XXX: Memory leak

		//solution->print();
		//ASSERT(solution->condition == NULL);

		CMatrix<Rational>* solM = solution->getBigParmSolution(0);
		//solM->print();
		if (solM == NULL) {
			delete solution;

			return NULL;
		} else {
			CMatrix<Rational> solC = solM->getSubMatrixColumns(solM->getColumnCount() - 1, 1);
			ASSERT(solC.getColumnCount() == 1);

			CVector<Rational>* result = new CVector<Rational>(solC.getRowCount());

			for (int q = 0; q < solM->getRowCount(); q++) (*result)[q] = solC[q][0];

			delete solution;

			return result;
		}
	}

	template <class R>
	CVector<typename R::RationalType>* pipVectorToCVector(PipVector * pv) {
		if (pv == NULL) return NULL;

		typedef typename R::RationalType Rational;
		CVector<Rational>& result = *(new CVector<Rational>(pv->nb_elements));

		for (int q = 0; q < pv->nb_elements; q++) {
			result[q] = Rational(CInteger(pv->the_vector[q]), CInteger(pv->the_deno[q]));
		}

		return &result;
	}

	template <class R>
	CMatrix<typename R::RationalType>* pipListToCMatrix(PipList * pl) {
		if (pl == NULL) return NULL;

		typedef typename R::RationalType Rational;
		vector<CVector<Rational>* > rows;

		while (pl != NULL) {
			rows.push_back(pipVectorToCVector<R>(pl->vector));
			pl = pl->next;
		}

		CMatrix<Rational>* result = new CMatrix<Rational>(rows.size() > 0 ? rows[0]->getLength() : -1);
		for (unsigned int q = 0; q < rows.size(); q++) {
			result->addRow(*rows[q]);
			delete rows[q];
		}

		return result;
	}

	template <class R>
	CPipQuast<R>* CPipQuast<R>::fromPipQuast(PipQuast* iOriginal) {
		if (iOriginal != NULL) {
			CPipQuast<R>& result = *(new CPipQuast<R>());

			result.nonNegativeCondQuast = CPipQuast<R>::fromPipQuast(iOriginal->next_then);
			result.negativeCondQuast = CPipQuast<R>::fromPipQuast(iOriginal->next_else);
			result.condition = pipVectorToCVector<R>(iOriginal->condition);
			result.solution = pipListToCMatrix<R>(iOriginal->list);

			return &result;
		} else return NULL;
	}
}

#endif
