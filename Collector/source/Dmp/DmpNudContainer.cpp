#include "Dmp/DmpNudContainer.h"

void DmpNudContainer::scanNudHits(const std::shared_ptr<DmpEvtNudRaw> &nudraw)
{
	for (int nud_ch = 0; nud_ch < DAMPE_NUD_channels; ++nud_ch)
	{
		adc[nud_ch] = (nudraw->fADC)[nud_ch];
		total_adc += adc[nud_ch];
		if (!nud_ch)
		{
			max_adc = adc[nud_ch];
			max_channel_id = nud_ch;
		}
		else
		{
			if (adc[nud_ch] > max_adc)
			{
				max_adc = adc[nud_ch];
				max_channel_id = nud_ch;
			}
		}
	}
}

const std::vector<int> DmpNudContainer::GetADC()
{
	return adc;
}
const int DmpNudContainer::GetTotalADC()
{
	return total_adc;
}
const int DmpNudContainer::GetMaxADC()
{
	return max_adc;
}
const int DmpNudContainer::GetMaxChannelID()
{
	return max_channel_id;
}

void DmpNudContainer::Reset()
{
	adc 					= std::vector<int>(DAMPE_NUD_channels, -999);
	total_adc 				= 0;
	max_adc 				= -999;
	max_channel_id 			= -999;
}