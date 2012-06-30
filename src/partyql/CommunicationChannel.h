#ifndef COMMUNICATIONCHANNEL_H_
#define COMMUNICATIONCHANNEL_H_

namespace AlgoTrans {
	template <class R> class CComputation;
	template <class R> class CStatement;
	template <class R> class CReference;
	template <class R> class CVariable;
	
	template <class R> class CModule;
	template <class R> class CAffineTransformation;
	template <class R> class CLattice;

	template <class R> class CCommunicationChannel {
		typedef CCommunicationChannel<R> CommunicationChannel;
		typedef CReference<R> Reference;
	private:
		const Reference& producingReference;
		const Reference& consumingReference;
		
		CPolyhedralDomain<R>* communicationSpace;
		
		CCommunicationChannel();
		
	public:
		CCommunicationChannel(const Reference& iProducingReference, const Reference& iConsumingReference);
		
		const CPolyhedralDomain<R>& getCommunicationSpace() const { return *communicationSpace; }
	};
}

#include "CommunicationChannel.cpp"

#endif /*COMMUNICATIONCHANNEL_H_*/
