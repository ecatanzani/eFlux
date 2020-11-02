#ifndef DMPNUDCONTAINER_H
#define DMPNUDCONTAINER_H

#include <vector>
#include <memory>

#include "DmpEvtNudRaw.h"
#include "DAMPE_geo_structure.h"

class DmpNudContainer
{
public:
	DmpNudContainer() : adc(DAMPE_NUD_channels, -999){};
	~DmpNudContainer(){};
	void scanNudHits(const std::shared_ptr<DmpEvtNudRaw> &nudraw);
	const std::vector<double> GetADC();
	const double GetTotalADC();
	const double GetMaxADC();
	const int GetMaxChannelID();
private:
	std::vector<double> adc;
	double total_adc = 0;
	double max_adc = -999;
	int max_channel_id = -999;
};

#endif