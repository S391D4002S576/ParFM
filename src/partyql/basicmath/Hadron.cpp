#ifdef HADRON_H_

namespace AlgoTrans {
	template <class DescriptorSet>
	CHadron<DescriptorSet>::CHadron(int iDimension, Description description)
	: constraints(NULL), generators(NULL) {
		dimension = iDimension;

		getDescriptorsP(description) = new DescriptorSet(dimension);
	};

	template <class DescriptorSet>
	CHadron<DescriptorSet>::CHadron(Description description, const DescriptorSet& descriptors)
	: constraints(NULL), generators(NULL) {
		dimension = descriptors.getDimension();

		getDescriptorsP(description) = new DescriptorSet(descriptors);
	};

	template <class DescriptorSet>
	CHadron<DescriptorSet>::CHadron(const CHadron& iOriginal)
	: constraints(NULL), generators(NULL) {
		*this = iOriginal;
	};

	template <class DescriptorSet>
	void CHadron<DescriptorSet>::cleanDescriptors(Description description) {
		DescriptorSetP& descriptors = getDescriptorsP(description);

		delete descriptors;
		descriptors = NULL;
	}

	template <class DescriptorSet>
	void CHadron<DescriptorSet>::cleanup() {
		cleanDescriptors(G);
		cleanDescriptors(C);
	}

	template <class DescriptorSet>
	CHadron<DescriptorSet>::~CHadron() {
		cleanup();
	}

	template <class DescriptorSet>
	CHadron<DescriptorSet>& CHadron<DescriptorSet>::operator = (const CHadron<DescriptorSet>& other) {
		if (this != &other) {
			cleanup();
			dimension = other.dimension;
			generators = (other.generators == NULL) ? NULL : new DescriptorSet(*other.generators);
			constraints = (other.constraints == NULL) ? NULL : new DescriptorSet(*other.constraints);
		}

		return *this;
	}

	template <class DescriptorSet>
	void CHadron<DescriptorSet>::setDescriptors(Description description, const DescriptorSet& descriptors) {
		ASSERT(dimension == descriptors.getDimension());

		if (getDescriptorsP(description) != &descriptors) {
			cleanup();
			getDescriptorsP(description) = new DescriptorSet(descriptors);
		}
	}

	template <class DescriptorSet>
	const DescriptorSet& CHadron<DescriptorSet>::getDescriptors(Description description) {
		if (getDescriptorsP(description) == NULL) { // Create it if we don't have it yet
			getDescriptorsP(description) = new DescriptorSet(!(*getDescriptorsP(!description)));
		}

		return *getDescriptorsP(description);
	}

	template <class DescriptorSet>
	typename DescriptorSet::Pointer& CHadron<DescriptorSet>::getDescriptorsP(Description description) {
		return (description == G) ? generators : constraints;
	}

	template <class DescriptorSet>
	string CHadron<DescriptorSet>::toString(bool forceBothDescriptions) {
		std::string result = getTypeName();
		if (forceBothDescriptions){
		    const DescriptorSet& gen = getDescriptors(G);
		    const DescriptorSet& con = getDescriptors(C);
			return result + "[G: " + gen.toString() + " | C: " + con.toString() + "]";
		} else {
			return result + "[G: " + ((generators != NULL) ? generators->toString() : "NULL")
		       + " | C: " + ((constraints != NULL) ? constraints->toString() : "NULL") + "]";
		}
	}

	template <class DescriptorSet>
	string CHadron<DescriptorSet>::toStringHtml(bool forceBothDescriptions) {
		std::string result = "<table><tr><td colspan=\"2\">" + getTypeNameHtml() + "</td></tr>";
		if (forceBothDescriptions || true){
		   const DescriptorSet& gen = getDescriptors(G);
		   const DescriptorSet& con = getDescriptors(C);
			return result + "<tr><td>G:</td><td>" + gen.toStringHtml()
			        + "</td></tr><tr><td>C:</td><td>" + con.toStringHtml() + "</td></tr></table>";
		} else {
			return result + "<tr><td>G:</td><td>" + ((generators != NULL) ? generators->toStringHtml() : "NULL")
			        + "</td></tr><tr><td>C:</td><td>" + ((constraints != NULL) ? constraints->toStringHtml() : "NULL") + "</td></tr></table>";
		}
	}

	template <class DescriptorSet>
	void CHadron<DescriptorSet>::dualize() {
		DescriptorSet* qt = generators;
		generators = constraints;
		constraints = qt;
	}

	template <class DescriptorSet>
	CHadron<DescriptorSet> CHadron<DescriptorSet>::getDual() const {
		CHadron<DescriptorSet> result = *this;

		result.dualize();

		return result;
	}
}

#endif
