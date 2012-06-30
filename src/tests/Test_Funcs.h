#ifndef TEST_FUNCS_H_
#define TEST_FUNCS_H_

#include <stdarg.h>

#include "../partyql/basicmath/scalar/Integer.h"
#include "../partyql/basicmath/Vektor.h"
#include "../partyql/basicmath/AffineTransformation.h"
#include "../partyql/basicmath/Polyheder.h"
#include "../partyql/basicmath/PolyhedralDomain.h"
#include "../partyql/basicmath/SetRelation.h"
#include "../partyql/basicmath/Lattice.h"
#include "../partyql/basicmath/Matrix.h"
#include "../partyql/basicmath/Module.h"
#include "../partyql/basicmath/Cone.h"

#include "../partyql/core/HadronDescription.h"
#include "../partyql/core/Descriptor.h"
#include "../partyql/core/Hadron.h"
#include "../partyql/core/Domain.h"
#include "../partyql/Computation.h"
#include "../partyql/HierarPart.h"
#include "../partyql/core/SetRelation.h"
#include "../partyql/core/util/Util.h"

using AlgoTrans::C;
using AlgoTrans::G;
using AlgoTrans::CVector;
using AlgoTrans::CInteger;
using AlgoTrans::CMatrix;
using AlgoTrans::CModule;
using AlgoTrans::CCone;
using AlgoTrans::CAffineTransformation;
using AlgoTrans::CModuleRelation;
using AlgoTrans::CLatticeRelation;
using AlgoTrans::CLattice;
using AlgoTrans::CPolyheder;
using AlgoTrans::CPolyhedralDomain;
using AlgoTrans::CFlat;

using AlgoTrans::Flat;
using AlgoTrans::Hadron;
using AlgoTrans::Description;
using AlgoTrans::DeskriptorSet;
using AlgoTrans::Deskriptor;
using AlgoTrans::CDisjunctiveDomain;
using AlgoTrans::SetRelation;

typedef CInteger I;
typedef CVector<I> IVector;
typedef CCone<I> ICone;
typedef CVector<int> IntVector;
typedef CMatrix<I> IMatrix;
typedef CModule<I> IModule;
typedef CFlat<IModule> IFlatModule;
typedef CAffineTransformation<I> IAffineTransformation;
typedef CModuleRelation<I>::T IModuleRelation;
typedef CLatticeRelation<I>::T ILatticeRelation;
typedef CLattice<I> ILattice;
typedef CPolyheder<I> IPolyheder;
typedef CPolyhedralDomain<I> IPolyhedralDomain;

typedef Deskriptor<I> IDeskriptor;
typedef DeskriptorSet<Deskriptor<I> > IDeskriptorSet;
typedef Hadron<DeskriptorSet<Deskriptor<I> > > IHadron;
typedef Flat<IHadron> Polyheder;
typedef SetRelation<Polyheder> ISetRelation;
typedef CDisjunctiveDomain<Polyheder> PolyhedralDomain;

struct IVector_POD { IVector* v; };
struct IPolyheder_POD { IPolyheder* p; };

typedef AlgoTrans::CComputation<I> IComputation;
typedef AlgoTrans::CHierarFinitePart<IComputation, I> IHierarFinitePart;
typedef AlgoTrans::CHierarSetPart<IComputation, I> IHierarSetPart;
typedef AlgoTrans::CHierarScatterPart<IComputation, I> IHierarScatterPart;
typedef AlgoTrans::CHierarPart<IComputation, I> IHierarPart;

// Parms: ints
IVector IV(int argCount, ...);
IVector_POD IV_(int argCount, ...);
IVector_POD IV_(IVector& v);

// Parms: ints
IntVector IntV(int argCount, ...);

// Parms: IVector_POD
IAffineTransformation IAT(int argCount, ...);
IMatrix IMx(int argCount, ...);
IModule IMod_Con(int argCount, ...);
IFlatModule IFMod_Con(int affineness, int argCount, ...);
IModule IMod_Gen(int argCount, ...);
ICone ICone_Con(int argCount, ...);
ICone ICone_Gen(int argCount, ...);
IModule IMod(int argCount, ...); // Default = Generators
//IPolyheder IP(int argCount, ...);
IPolyhedralDomain* pIP_D(int argCount, ...);

PolyhedralDomain* pIPDC(bool bidir, int argCount, ...);
PolyhedralDomain IPDC(bool bidir, int argCount, ...);
Polyheder IPC(bool bidir, int argCount, ...);
Polyheder IP(Description d, int argCount, ...);
IHadron IHC(bool bidir, int argCount, ...);
IHadron IH(Description d, bool bidir, int argCount, ...);
ISetRelation ISR(bool bidir, Description desc, int constDim, int pDim, int qDim, int argCount, ...);
IDeskriptorSet IDS(bool bidir, int argCount, ...);

//PolyhedralDomain* pIPDG(int argCount, ...);
IPolyheder_POD IP_(int argCount, ...);
IPolyhedralDomain IPD(int argCount, ...);
IModuleRelation IMR(int pDim, int qDim, int constDim, int argCount, ...);
IModuleRelation IMR_Con(int pDim, int qDim, int constDim, int argCount, ...);
ILatticeRelation ILR(int pDim, int qDim, int constDim, int argCount, ...);

#endif /*TEST_FUNCS_H_*/
