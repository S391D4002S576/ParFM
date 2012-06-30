#include "Util.h"

std::vector<int> getDimSetComplement(std::vector<int> dims, int dimension) {
	std::vector<bool> complDim = std::vector<bool>(dimension, true);
	for (unsigned int q = 0; q < dims.size(); q++) complDim[dims[q]] = false;
	std::vector<int> resultDims;
	for (unsigned int q = 0; q < complDim.size(); q++) if (complDim[q]) resultDims.push_back(q);

	return resultDims;
}

std::vector<int> getInterval(int start, int length) {
	std::vector<int> resultDims;

	for (int q = 0; q < length; q++) resultDims.push_back(start + q);

	return resultDims;
}
