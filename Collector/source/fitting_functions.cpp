#include "acceptance_fit.h"

#include "TMath.h"

double logisticFunction_0(double *x, double *par)
{
    double xx = log10(x[0]);
    double logistic = par[0] / (1 + TMath::Exp(-par[1]*xx));
    return logistic;
}

double logisticFunction_1(double *x, double *par)
{
    double xx = log10(x[0]);
    double logistic = par[0] / (1 + TMath::Exp(-par[1] * xx + par[2]*TMath::Power(xx,2) + par[3]));
    return logistic;
}

double logisticFunction_2(double *x, double *par)
{
    double xx = log10(x[0]);
    double logistic = (par[0] + par[4] * xx + par[5] * TMath::Power(xx,2) ) / (1 + TMath::Exp(-par[1] * xx + par[2]*TMath::Power(xx,2) + par[3]));
    return logistic;
}

double logisticFunction_3(double *x, double *par)
{
    double xx = log10(x[0]);
    double logistic = (par[0] + par[4] * xx + par[5] * TMath::Power(xx,2)) / (1 + TMath::Exp(-par[1] * xx + par[2]*TMath::Power(xx,2) + par[3] + par[6]*TMath::Power(xx,3) + par[7]*TMath::Power(xx,4) ));
    return logistic;
}

double logisticFunction_4(double *x, double *par)
{
    double xx = log10(x[0]);
    double logistic = (par[0] + par[4] * xx + par[5] * TMath::Power(xx,2) + par[8]*TMath::Power(xx,3) + par[9]*TMath::Power(xx,4)) / (1 + TMath::Exp(-par[1] * xx + par[2]*TMath::Power(xx,2) + par[3] + par[6]*TMath::Power(xx,3) + par[7]*TMath::Power(xx,4) ));
    return logistic;
}

double logisticFunction_5(double *x, double *par)
{
    double xx = log10(x[0]);
    double logistic = (par[0] + par[4] * xx + par[5] * TMath::Power(xx,2) + par[10]*TMath::Power(xx,4)) / (1 + TMath::Exp(-par[1] * xx + par[2]*TMath::Power(xx,2) + par[3] + par[6]*TMath::Power(xx,3) + par[7]*TMath::Power(xx,4) + par[8]*TMath::Power(xx,6) + par[9]*TMath::Power(xx,10)));
    return logistic;
}

double logisticFunction_6(double *x, double *par)
{
    double xx = log10(x[0]);
    double logistic = (par[0] + par[4] * xx + par[5] * TMath::Power(xx,2) + par[10]*TMath::Power(xx,4) + par[11]*TMath::Power(xx,6) + par[12]*TMath::Power(xx,8)) / (1 + TMath::Exp(-par[1] * xx + par[2]*TMath::Power(xx,2) + par[3] + par[6]*TMath::Power(xx,3) + par[7]*TMath::Power(xx,4) + par[8]*TMath::Power(xx,6) + par[9]*TMath::Power(xx,10)));
    return logistic;
}

double logisticFunction_7(double *x, double *par)
{
    double xx = log10(x[0]);
    double logistic = (par[0] + par[4] * xx + par[5] * TMath::Power(xx,2) + par[10]*TMath::Power(xx,4) + par[11]*TMath::Power(xx,6) + par[12]*TMath::Power(xx,8) + par[13]*TMath::Power(xx,10) + par[14]*TMath::Power(xx,12)) / (1 + TMath::Exp(-par[1] * xx + par[2]*TMath::Power(xx,2) + par[3] + par[6]*TMath::Power(xx,3) + par[7]*TMath::Power(xx,4) + par[8]*TMath::Power(xx,6) + par[9]*TMath::Power(xx,10)));
    return logistic;
}

