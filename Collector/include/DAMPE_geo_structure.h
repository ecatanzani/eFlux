#ifndef DAMPE_GEO_STRUCTURE_H
#define DAMPE_GEO_STRUCTURE_H

// DAMPE struct
const int DAMPE_psd_nLayers = 2;
const int DAMPE_bgo_nLayers = 14;
const int nSTKladders = 192;
const double BGO_TopZ = 46;
const double BGO_BottomZ = 448;
const double BGO_SideXY = 301.25;
const double BGO_layer_displacement = (BGO_BottomZ-BGO_TopZ)/(DAMPE_bgo_nLayers-1);

#endif