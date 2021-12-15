#ifndef BDTINFO_H
#define BDTINFO_H

#include "main.h"

struct bdt_vars {
    float rmslayer_norm_1 {0};
    float rmslayer_norm_2 {0};
    float rmslayer_norm_3 {0};
    float rmslayer_norm_4 {0};
    float rmslayer_norm_5 {0};
    float rmslayer_norm_6 {0};
    float rmslayer_norm_7 {0};
    float rmslayer_norm_8 {0};
    float rmslayer_norm_9 {0};
    float rmslayer_norm_10 {0};
    float rmslayer_norm_11 {0};
    float rmslayer_norm_12 {0};
    float rmslayer_norm_13 {0};
    float rmslayer_norm_14 {0};
    float fraclayer_norm_1 {0};
    float fraclayer_norm_2 {0};
    float fraclayer_norm_3 {0};
    float fraclayer_norm_4 {0};
    float fraclayer_norm_5 {0};
    float fraclayer_norm_6 {0};
    float fraclayer_norm_7 {0};
    float fraclayer_norm_8 {0};
    float fraclayer_norm_9 {0};
    float fraclayer_norm_10 {0};
    float fraclayer_norm_11 {0};
    float fraclayer_norm_12 {0};
    float fraclayer_norm_13 {0};
    float fraclayer_norm_14 {0};
    float sumrms_norm {0};
    float fraclastlayer_norm {0};
    float xtrl_norm {0};
    float xtrl {0};
};

struct data_vars {
    double rmslayer_norm_1 {0};
    double rmslayer_norm_2 {0};
    double rmslayer_norm_3 {0};
    double rmslayer_norm_4 {0};
    double rmslayer_norm_5 {0};
    double rmslayer_norm_6 {0};
    double rmslayer_norm_7 {0};
    double rmslayer_norm_8 {0};
    double rmslayer_norm_9 {0};
    double rmslayer_norm_10 {0};
    double rmslayer_norm_11 {0};
    double rmslayer_norm_12 {0};
    double rmslayer_norm_13 {0};
    double rmslayer_norm_14 {0};
    double fraclayer_norm_1 {0};
    double fraclayer_norm_2 {0};
    double fraclayer_norm_3 {0};
    double fraclayer_norm_4 {0};
    double fraclayer_norm_5 {0};
    double fraclayer_norm_6 {0};
    double fraclayer_norm_7 {0};
    double fraclayer_norm_8 {0};
    double fraclayer_norm_9 {0};
    double fraclayer_norm_10 {0};
    double fraclayer_norm_11 {0};
    double fraclayer_norm_12 {0};
    double fraclayer_norm_13 {0};
    double fraclayer_norm_14 {0};
    double sumrms_norm {0};
    double fraclastlayer_norm {0};
    double xtrl_norm {0};
    double xtrl {0};

    double evt_corr_energy {0};
};

extern void ExtractBDTInfo(in_args input_args);

#endif