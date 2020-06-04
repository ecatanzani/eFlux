#include "acceptance_fit.h"

#include "TMath.h"

double logisticFunction_0(double *x, double *par)
{
    double xx = log10(x[0] + 1);
    double logistic = (par[0]) / (1 + TMath::Exp(-par[1] * (par[2] * xx + par[3])));
    return logistic;
}

double logisticFunction_1(double *x, double *par)
{
    double xx = log10(x[0] + 1);
    double logistic = (par[0] + par[4] * xx + par[5] * TMath::Power(xx, 2)) / (1 + TMath::Exp(-par[1] * (par[2] * xx + par[3])));
    return logistic;
}

double logisticFunction_2(double *x, double *par)
{
    double xx = log10(x[0] + 1);
    double logistic = (par[0] + par[4] * xx + par[5] * TMath::Power(xx, 2)) / (1 + TMath::Exp(-par[1] * (par[2] * xx + par[3] * par[6] * TMath::Power(xx, 2))));
    return logistic;
}

double logisticFunction_3(double *x, double *par)
{
    double xx = log10(x[0] + 1);
    double logistic = (par[0] + par[4] * xx + par[5] * TMath::Power(xx, 2)) / (1 + TMath::Exp(-par[1] * (par[2] * xx + par[3] * par[6] * TMath::Power(xx, 2) + par[7] * TMath::Power(xx, par[8]))));
    return logistic;
}

double logisticFunction_4(double *x, double *par)
{
    double xx = log10(x[0] + 1);
    double logistic = (par[0] + par[4] * xx + par[5] * TMath::Power(xx, 2) + par[9] * TMath::Power(xx, par[10])) / (1 + TMath::Exp(-par[1] * (par[2] * xx + par[3] * par[6] * TMath::Power(xx, 2) + par[7] * TMath::Power(xx, par[8]))));
    return logistic;
}

double logisticFunction_5(double *x, double *par)
{
    double xx = log10(x[0] + 1);
    double logistic = (par[0] + par[4] * xx + par[5] * TMath::Power(xx, 2) + par[11] * TMath::Power(xx, 3) + par[12] * TMath::Power(xx, 4) + par[9] * TMath::Power(xx, par[10])) / (1 + TMath::Exp(-par[1] * (par[2] * xx + par[3] * par[6] * TMath::Power(xx, 2) + par[7] * TMath::Power(xx, par[8]))));
    return logistic;
}

double logisticFunction_6(double *x, double *par)
{
    double xx = log10(x[0] + 1);
    double logistic = (par[0] + par[4] * xx + par[5] * TMath::Power(xx, 2) + par[11] * TMath::Power(xx, 3) + par[12] * TMath::Power(xx, 4) + par[13] * TMath::Power(xx, 5) + par[14] * TMath::Power(xx, 6) + par[9] * TMath::Power(xx, par[10])) / (1 + TMath::Exp(-par[1] * (par[2] * xx + par[3] * par[6] * TMath::Power(xx, 2) + par[7] * TMath::Power(xx, par[8]))));
    return logistic;
}

double logisticFunction_7(double *x, double *par)
{
    double xx = log10(x[0] + 1);
    double logistic = (par[0] + par[4] * xx + par[5] * TMath::Power(xx, 2) + par[11] * TMath::Power(xx, 3) + par[12] * TMath::Power(xx, 4) + par[13] * TMath::Power(xx, 5) + par[14] * TMath::Power(xx, 6) + par[15] * TMath::Power(xx, 7) + par[16] * TMath::Power(xx, 8) + par[17] * TMath::Power(xx, 9) + par[9] * TMath::Power(xx, par[10])) / (1 + TMath::Exp(-par[1] * (par[2] * xx + par[3] * par[6] * TMath::Power(xx, 2) + par[7] * TMath::Power(xx, par[8]))));
    return logistic;
}

double logisticFunction_8(double *x, double *par)
{
    double xx = log10(x[0] + 1);
    double logistic = (par[0] + par[4] * xx + par[5] * TMath::Power(xx, 2) + par[11] * TMath::Power(xx, 3) + par[12] * TMath::Power(xx, 4) + par[13] * TMath::Power(xx, 5) + par[14] * TMath::Power(xx, 6) + par[15] * TMath::Power(xx, 7) + par[16] * TMath::Power(xx, 8) + par[17] * TMath::Power(xx, 9) + par[18] * TMath::Power(xx, 10) + par[9] * TMath::Power(xx, par[10])) / (1 + TMath::Exp(-par[1] * (par[2] * xx + par[3] * par[6] * TMath::Power(xx, 2) + par[19] * TMath::Power(xx, 3) + par[20] * TMath::Power(xx, 4) + par[7] * TMath::Power(xx, par[8]))));
    return logistic;
}
