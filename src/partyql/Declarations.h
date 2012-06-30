#ifndef DECLARATIONS_H_
#define DECLARATIONS_H_

#define ITTT(IteratorType, IteratorVariable, IteratorSet) for (IteratorType IteratorVariable = IteratorSet.begin(); IteratorVariable != IteratorSet.end(); ++IteratorVariable)
#define ITT(IteratorSetType, IteratorVariable, IteratorSet) ITTT(typename IteratorSetType::iterator, IteratorVariable, (IteratorSet))
#define ITTc(IteratorSetType, IteratorVariable, IteratorSet) ITTT(typename IteratorSetType::const_iterator, IteratorVariable, (IteratorSet))

#endif
