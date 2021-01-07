#include "xtrX_computation.h"

double xtrX_computation(
	const double sumRms,
	const double lastFracLayer)
{
	if (lastFracLayer != -1)
		return 0.125e-6 * pow(sumRms, 4) * lastFracLayer;
	else
		return -999;
}