#ifndef DATA_TMPSTRUCT_H
#define DATA_TMPSTRUCT_H

#include <memory>

#include "energy.h"

struct _tmp_energy_data
{
	double raw;
	double correct;
	double correct_w;
};

extern std::shared_ptr<_tmp_energy_data> fillDataEnergyTmpStruct(energy &evt_energy);

#endif