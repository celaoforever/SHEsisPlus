/*
 * Haplotype.h
 *
 *  Created on: Aug 19, 2014
 *      Author: ada
 */

#ifndef HAPLOTYPE_H_
#define HAPLOTYPE_H_
#include "SHEsisData.h"
#include "utility.h"
#include "ArrayStorage.h"
#include "IndexingVariables.h"
namespace SHEsis {
class Haplotype {
public:
	Haplotype(SHEsisData data):data(data){};
	virtual ~Haplotype();
	SHEsisData& data;
};

} /* namespace SHEsis */

#endif /* HAPLOTYPE_H_ */
