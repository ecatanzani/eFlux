#include "DmpNudContainer.h"

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

const std::vector<double> DmpNudContainer::GetADC()
{
	return adc;
}
const double DmpNudContainer::GetTotalADC()
{
	return total_adc;
}
const double DmpNudContainer::GetMaxADC()
{
	return max_adc;
}
const int DmpNudContainer::GetMaxChannelID()
{
	return max_channel_id;
}