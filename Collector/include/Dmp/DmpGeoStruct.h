#ifndef DMPGEOSTRUCT_H
#define DMPGEOSTRUCT_H

// DAMPE struct
const int DAMPE_psd_nLayers             {2};
const int DAMPE_bgo_nLayers             {14};
const int DAMPE_bgo_bars_layer          {22};
const int DAMPE_stk_planes              {6};
const int DAMPE_NUD_channels            {4};
const int nSTKladders                   {192};
const double BGO_TopZ                   {46};
const double BGO_BottomZ                {448};
const double BGO_SideXY                 {301.25};
const double BGO_layer_displacement     {(BGO_BottomZ-BGO_TopZ)/(DAMPE_bgo_nLayers-1)};
const double PSD_TopZ                   {-331.7};
const double PSD_SideXY                 {410};
const double STK_TopZ                   {-210};
const double STK_BottomZ                {-41.968};
const double STK_SideXY                 {400};

#endif