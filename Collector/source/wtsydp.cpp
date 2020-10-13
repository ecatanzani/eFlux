#include "wtsydp.h"

double wtsydp(
    const float minene,
    const float maxene,
    const float index)
{
    float dene = maxene - minene;
    if (index != -1)
        return pow(fabs((pow(maxene, index + 1) - pow(minene, index + 1)) / ((index + 1) * dene)), 1. / index);
    else
        return dene / log(maxene / minene);
}
