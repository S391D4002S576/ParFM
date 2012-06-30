#include "Test_Funcs.h"

#include "../partyql/basicmath/Flat.h"
#include "../partyql/basicmath/Hadron.h"
#include "../partyql/basicmath/VectorDescriptor.h"
#include "../partyql/core/Vektor.h"
using AlgoTrans::Vektor;

IVector_POD IV_(IVector& v) {
	IVector_POD pod;

	pod.v = new IVector(v);

	return pod;
};

IVector IV(int argCount, ...) {
	va_list argList;
	va_start(argList, argCount);

	IVector result = IVector();
	while (argCount-- > 0) result.appendElement(va_arg(argList, int));

	va_end(argList);

	return result;
}

IVector_POD IV_(int argCount, ...) {
	va_list argList;
	va_start(argList, argCount);

	IVector_POD result;
	result.v = new IVector();
	while (argCount-- > 0) result.v->appendElement(va_arg(argList, int));

	va_end(argList);

	return result;
}

IntVector IntV(int argCount, ...) {
	va_list argList;
	va_start(argList, argCount);

	IntVector result = IntVector();
	while (argCount-- > 0) result.appendElement(va_arg(argList, int));

	va_end(argList);

	return result;
}

IMatrix extractMatrix(const vector<IVector*>& vs) {
	IMatrix result = (vs.size() > 0) ? IMatrix(vs[0]->getLength()) : IMatrix(-1);

	for (unsigned int q = 0; q < vs.size(); q++) {
		result.addRow(*vs[q]);
		ASSERT(vs[q]->getLength() == result.getColumnCount());
	}

	return result;
}

IPolyhedralDomain extractPolyhedralDomain(const vector<IPolyheder*>& vs) {
	IPolyhedralDomain result = (vs.size() > 0) ? IPolyhedralDomain(vs[0]->getSpaceDimension()) : IPolyhedralDomain(-1);

	for (unsigned int q = 0; q < vs.size(); q++) {
		result.addPolyheder(*vs[q]);
		ASSERT(vs[q]->getSpaceDimension() == result.getSpaceDimension());
	}

	return result;
}

IMatrix IMx(int argCount, ...) {
	va_list argList;
	va_start(argList, argCount);

	vector<IVector*> vs = vector<IVector*>();
	while (argCount-- > 0) vs.push_back(va_arg(argList, IVector_POD).v);
	va_end(argList);

	return IMatrix(extractMatrix(vs));
}

IAffineTransformation IAT(int argCount, ...) {
	va_list argList;
	va_start(argList, argCount);

	vector<IVector*> vs = vector<IVector*>();
	while (argCount-- > 0) vs.push_back(va_arg(argList, IVector_POD).v);
	va_end(argList);

	return IAffineTransformation(extractMatrix(vs));
}

IModule IMod_Con(int argCount, ...) {
	va_list argList;
	va_start(argList, argCount);

	if (argCount == 0) {
		int dim = va_arg(argList, int);
		va_end(argList);
		return IModule(dim, C);
	} else {
		vector<IVector*> vs = vector<IVector*>();
		while (argCount-- > 0) vs.push_back(va_arg(argList, IVector_POD).v);
		va_end(argList);
		return IModule(C, matrixToDescriptorSet(extractMatrix(vs)));
	}
}

IFlatModule IFMod_Con(int affineness, int argCount, ...) {
	va_list argList;
	va_start(argList, argCount);

	if (argCount == 0) {
		int dim = va_arg(argList, int);
		va_end(argList);
		return IFlatModule(affineness, IModule(dim, C));
	} else {
		vector<IVector*> vs = vector<IVector*>();
		while (argCount-- > 0) vs.push_back(va_arg(argList, IVector_POD).v);
		va_end(argList);
		return IFlatModule(affineness, IModule(C, matrixToDescriptorSet(extractMatrix(vs))));
	}
}


IModule IMod_Gen(int argCount, ...) {
	va_list argList;
	va_start(argList, argCount);

	vector<IVector*> vs = vector<IVector*>();
	while (argCount-- > 0) vs.push_back(va_arg(argList, IVector_POD).v);
	va_end(argList);

	return IModule(G, matrixToDescriptorSet(extractMatrix(vs)));
}

ICone ICone_Con(int argCount, ...) {
	va_list argList;
	va_start(argList, argCount);

	if (argCount == 0) {
		int dim = va_arg(argList, int);
		va_end(argList);
		return ICone(dim, C);
	} else {
		vector<IVector*> vs = vector<IVector*>();
		while (argCount-- > 0) vs.push_back(va_arg(argList, IVector_POD).v);
		va_end(argList);
		return ICone(C, matrixToDescriptorSet(extractMatrix(vs)));
	}
}

ICone ICone_Gen(int argCount, ...) {
	va_list argList;
	va_start(argList, argCount);

	vector<IVector*> vs = vector<IVector*>();
	while (argCount-- > 0) vs.push_back(va_arg(argList, IVector_POD).v);
	va_end(argList);

	return ICone(G, matrixToDescriptorSet(extractMatrix(vs)));
}

IModule IMod(int argCount, ...) {
	va_list argList;
	va_start(argList, argCount);

	vector<IVector*> vs = vector<IVector*>();
	while (argCount-- > 0) vs.push_back(va_arg(argList, IVector_POD).v);
	va_end(argList);

	return IModule(G, matrixToDescriptorSet(extractMatrix(vs)));
}

IPolyheder_POD IP_(int argCount, ...) {
	va_list argList;
	va_start(argList, argCount);

	vector<IVector*> vs = vector<IVector*>();
	while (argCount-- > 0) vs.push_back(va_arg(argList, IVector_POD).v);
	va_end(argList);

	IPolyheder_POD result;
	result.p = new IPolyheder(extractMatrix(vs));
	return result;
}

IPolyheder IP(int argCount, ...) {
	va_list argList;
	va_start(argList, argCount);

	vector<IVector*> vs = vector<IVector*>();

	if (argCount == 0) {
		int dim = va_arg(argList, int);
		va_end(argList);
		return IPolyheder(dim);
	} else {
		while (argCount-- > 0) vs.push_back(va_arg(argList, IVector_POD).v);
		va_end(argList);
		return IPolyheder(extractMatrix(vs));
	}
}

IPolyhedralDomain* pIP_D(int argCount, ...) {
	va_list argList;
	va_start(argList, argCount);

	vector<IVector*> vs = vector<IVector*>();

	if (argCount == 0) {
		int dim = va_arg(argList, int);
		va_end(argList);
		return new IPolyhedralDomain(IPolyheder(dim));
	} else {
		while (argCount-- > 0) vs.push_back(va_arg(argList, IVector_POD).v);
		va_end(argList);
		return new IPolyhedralDomain(IPolyheder(extractMatrix(vs)));
	}
}

template <class R>
DeskriptorSet<Deskriptor<R> > matrixToDeskriptorSet(const CMatrix<R>& matrix, bool bidir) {
	DeskriptorSet<Deskriptor<R> > result = DeskriptorSet<Deskriptor<R> >(matrix.getColumnCount());

	Vektor<R> v = Vektor<R>(matrix.getColumnCount());

	for (int q = 0; q < matrix.getRowCount(); q++) {
		for (int r = 0; r < v.getLength(); r++) v[r] = matrix[q][r];
   		result.addDeskriptor(Deskriptor<R>(v, bidir));
	}

	return result;
}



PolyhedralDomain* pIPDC(int argCount, ...) {
	va_list argList;
	va_start(argList, argCount);

	vector<IVector*> vs = vector<IVector*>();

	if (argCount == 0) {
		int dim = va_arg(argList, int);
		va_end(argList);
		return new PolyhedralDomain(Polyheder::universe(dim));
	} else {
		while (argCount-- > 0) vs.push_back(va_arg(argList, IVector_POD).v);
		va_end(argList);
		return new PolyhedralDomain(Polyheder(C, matrixToDeskriptorSet(extractMatrix(vs), false)));
	}
}

PolyhedralDomain IPDC(bool bidir, int argCount, ...) {
	va_list argList;
	va_start(argList, argCount);

	vector<IVector*> vs = vector<IVector*>();

	if (argCount == 0) {
		int dim = va_arg(argList, int);
		va_end(argList);
		return PolyhedralDomain(Polyheder::universe(dim));
	} else {
		while (argCount-- > 0) vs.push_back(va_arg(argList, IVector_POD).v);
		va_end(argList);
		return PolyhedralDomain(Polyheder(C, matrixToDeskriptorSet(extractMatrix(vs), bidir)));
	}
}

Polyheder IPC(bool bidir, int argCount, ...) {
	va_list argList;
	va_start(argList, argCount);

	vector<IVector*> vs = vector<IVector*>();

	if (argCount == 0) {
		int dim = va_arg(argList, int);
		va_end(argList);
		return Polyheder::universe(dim);
	} else {
		while (argCount-- > 0) vs.push_back(va_arg(argList, IVector_POD).v);
		va_end(argList);
		return Polyheder(C, matrixToDeskriptorSet(extractMatrix(vs), false));
	}
}

Polyheder IP(Description d, int argCount, ...) {
	va_list argList;
	va_start(argList, argCount);

	vector<IVector*> vs = vector<IVector*>();

	if (argCount == 0) {
		int dim = va_arg(argList, int);
		va_end(argList);
		return (d == C) ? Polyheder::universe(dim) : Polyheder::source(dim);
	} else {
		while (argCount-- > 0) vs.push_back(va_arg(argList, IVector_POD).v);
		va_end(argList);
		return Polyheder(d, matrixToDeskriptorSet(extractMatrix(vs), false));
	}
}

IDeskriptorSet IDS(bool bidir, int argCount, ...) {
	va_list argList;
	va_start(argList, argCount);

	vector<IVector*> vs = vector<IVector*>();

	if (argCount == 0) {
		int dim = va_arg(argList, int);
		va_end(argList);
		return IDeskriptorSet(dim);
	} else {
		while (argCount-- > 0) vs.push_back(va_arg(argList, IVector_POD).v);
		va_end(argList);
		return matrixToDeskriptorSet(extractMatrix(vs), bidir);
	}
}

IHadron IHC(bool bidir, int argCount, ...) {
	va_list argList;
	va_start(argList, argCount);

	vector<IVector*> vs = vector<IVector*>();

	if (argCount == 0) {
		int dim = va_arg(argList, int);
		va_end(argList);
		return IHadron::universe(dim);
	} else {
		while (argCount-- > 0) vs.push_back(va_arg(argList, IVector_POD).v);
		va_end(argList);
		return IHadron(C, matrixToDeskriptorSet(extractMatrix(vs), bidir));
	}
}

IHadron IH(Description d, bool bidir, int argCount, ...) {
	va_list argList;
	va_start(argList, argCount);

	vector<IVector*> vs = vector<IVector*>();

	if (argCount == 0) {
		int dim = va_arg(argList, int);
		va_end(argList);
		return (d == C) ? IHadron::universe(dim) : IHadron::source(dim);
	} else {
		while (argCount-- > 0) vs.push_back(va_arg(argList, IVector_POD).v);
		va_end(argList);
		return IHadron(d, matrixToDeskriptorSet(extractMatrix(vs), bidir));
	}
}

ISetRelation ISR(bool bidir, Description desc, int constDim, int pDim, int qDim, int argCount, ...) {
	va_list argList;
	va_start(argList, argCount);

	vector<IVector*> vs = vector<IVector*>();

	if (argCount == 0) {
		return ISetRelation(constDim, pDim, qDim, (desc == C) ? Polyheder::universe(constDim + pDim + qDim) : Polyheder::source(constDim + pDim + qDim));
	} else {
		while (argCount-- > 0) vs.push_back(va_arg(argList, IVector_POD).v);
		va_end(argList);
		return ISetRelation(constDim, pDim, qDim, Polyheder(desc, matrixToDeskriptorSet(extractMatrix(vs), bidir)));
	}
}

IPolyhedralDomain IPD(int argCount, ...) {
	va_list argList;
	va_start(argList, argCount);

	vector<IPolyheder*> vs = vector<IPolyheder*>();
	while (argCount-- > 0) vs.push_back(va_arg(argList, IPolyheder_POD).p);
	va_end(argList);

	return extractPolyhedralDomain(vs);
}

ILatticeRelation ILR(int pDim, int qDim, int constDim, int argCount, ...) {
	va_list argList;
	va_start(argList, argCount);

	if (argCount == 0) {
		va_end(argList);
		return ILatticeRelation(pDim, qDim, CFlat<ILattice>(constDim, ILattice(pDim + qDim + constDim, G)));
	} else {
		vector<IVector*> vs = vector<IVector*>();
		while (argCount-- > 0) vs.push_back(va_arg(argList, IVector_POD).v);
		va_end(argList);

		return ILatticeRelation(pDim, qDim, CFlat<ILattice>(constDim, ILattice(G, matrixToDescriptorSet(extractMatrix(vs)))));
	}
}

IModuleRelation IMR(int pDim, int qDim, int constDim, int argCount, ...) {
	va_list argList;
	va_start(argList, argCount);

	if (argCount == 0) {
		va_end(argList);
		return IModuleRelation(pDim, qDim, CFlat<IModule>(constDim, IModule(pDim + qDim + constDim, G)));
	} else {
		vector<IVector*> vs = vector<IVector*>();
		while (argCount-- > 0) vs.push_back(va_arg(argList, IVector_POD).v);
		va_end(argList);

		return IModuleRelation(pDim, qDim, CFlat<IModule>(constDim, IModule(G, matrixToDescriptorSet(extractMatrix(vs)))));
	}
}

IModuleRelation IMR_Con(int pDim, int qDim, int constDim, int argCount, ...) {
	va_list argList;
	va_start(argList, argCount);

	if (argCount == 0) {
		va_end(argList);
		return IModuleRelation(pDim, qDim, CFlat<IModule>(constDim, IModule(pDim + qDim + constDim, C)));
	} else {
		vector<IVector*> vs = vector<IVector*>();
		while (argCount-- > 0) vs.push_back(va_arg(argList, IVector_POD).v);
		va_end(argList);
		return IModuleRelation(pDim, qDim, CFlat<IModule>(constDim, IModule(C, matrixToDescriptorSet(extractMatrix(vs)))));
	}
}

