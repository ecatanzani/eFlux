#ifndef DMPNUDCONTAINER_H
#define DMPNUDCONTAINER_H

#include <vector>
#include <memory>

#include "DmpEvtNudRaw.h"
#include "DAMPE_geo_structure.h"

class DmpNudContainer
{
public:
	DmpNudContainer() : adc(4, 0){};
	~DmpNudContainer(){};
	void scanNudHits(const std::shared_ptr<DmpEvtNudRaw> &nudraw);
	const std::vector<double> GetADC();
	const double GetTotalADC();
	const double GetMaxADC();
	const int GetMaxChannelID();
private:
	std::vector<double> adc;
	double total_adc = 0;
	double max_adc = 0;
	int max_channel_id = 0;
};

#endif