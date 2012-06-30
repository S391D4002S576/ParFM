#ifndef STRUCTURES_H_
#define STRUCTURES_H_

#include <vector>
#include <list>

#include "Declarations.h"

namespace AlgoTrans {
	/*
	 * Reqs:
	 * - ParentSetClass
	 *   - iterator, begin(), end()
	 *   - size()
	 * - ElementType
	 *   - getIndex()
	 * Use: Parent set is a fixed set of which we want to temporarily consider subsets
	 */
	template <class ParentSetClass, class ElementType> class CFastSubSet {
		typedef typename std::list<ElementType* >::iterator ElementListIterator;
	private:
		std::vector<ElementListIterator* > iterators;
		std::vector<bool> included;
		std::list<ElementType* > includeList, excludeList;
		const ParentSetClass& parentSet;

	public:
		bool isIncluded(int index) const { return included[index]; }
		bool isIncluded(const ElementType& element) const { return isIncluded(element.getIndex()); }

		CFastSubSet(const ParentSetClass& iParentSet, bool initiallyEmpty) : parentSet(iParentSet) {
			iterators = vector<ElementListIterator* >(parentSet.size(), NULL);
			int q = 0;
			ITTc(ParentSetClass, parentI, parentSet) {
				if (!initiallyEmpty) {
					includeList.push_front(*parentI);
					iterators[q] = new ElementListIterator(includeList.begin());
				} else {
					excludeList.push_front(*parentI);
					iterators[q] = new ElementListIterator(excludeList.begin());
				}
				included.push_back(!initiallyEmpty);
				q++;
			}
		}

		std::vector<ElementType* > includedsAsVector() {
			std::vector<ElementType* > result;

			ITT(std::list<ElementType* >, eli, includeList) result.push_back(*eli);

			return result;
		}

		std::vector<ElementType* > excludedsAsVector() {
			std::vector<ElementType* > result;

			ITT(std::list<ElementType* >, eli, excludeList) result.push_back(*eli);

			return result;
		}

		void includeElement(ElementType& element) { int elix = element.getIndex();
			if (!isIncluded(elix)) {
				excludeList.erase(*iterators[elix]);
				delete iterators[elix];

				includeList.push_front(&element);
				iterators[elix] = new ElementListIterator(includeList.begin());
			}
		}

		void excludeElement(ElementType& element) { int elix = element.getIndex();
			if (isIncluded(elix)) {
				includeList.erase(*iterators[elix]);
				delete iterators[elix];

				excludeList.push_front(&element);
				iterators[elix] = new ElementListIterator(excludeList.begin());
			}
		}

		~CFastSubSet() {
			ITT(std::vector<ElementListIterator* >, i, iterators) delete *i;
		}
	};

	template <class Key, class Value> class CIxMap {
	private:
		std::vector<Value* > elements;
	public:
		Value* operator [] (const Key& key) { return elements[key->getIndex()]; };

		CIxMap(int setSize) { elements = vector<Value* >(setSize, NULL); }

		~CIxMap() { ITT(std::vector<Value* >, e, elements) delete *e; }
	};

	template <class Key, class Value, class ValueSet> class CIxLink {
	private:
		ValueSet values;
	public:
		Value& operator [] (const Key& key) { return values(key->getIndex()); };

		CIxLink(ValueSet& iValues) : values(iValues) { }
	};
}

#endif
