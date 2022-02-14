#include "Data/data_tmpstruct.h"

std::shared_ptr<_tmp_energy_data> fillDataEnergyTmpStruct(energy &evt_energy)
{
    std::shared_ptr<_tmp_energy_data> _energy_res = std::make_shared<_tmp_energy_data>();

	_energy_res->raw = evt_energy.GetRawEnergy();
	_energy_res->correct = evt_energy.GetCorrEnergy();
	_energy_res->correct_w = evt_energy.GetCorrWeight();

	return _energy_res;
}