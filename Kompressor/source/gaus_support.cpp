#include "gaussianize.h"

#include "TF1.h"

void extract_lamda_info(
    const std::vector<std::vector<std::vector<std::shared_ptr<TH1D>>>> input_histo,
    std::vector<double> &goodness,
    std::vector<double> &best_lambda,
    std::vector<int> &best_lambda_idx,
    const double lambda_start,
    const double lambda_step,
    const int bin_idx,
    const int lambda_idx,
    const int layer)
{
    if (input_histo[bin_idx-1][lambda_idx][layer]->GetEntries()>5)
    {
        auto xmin = input_histo[bin_idx-1][lambda_idx][layer]->GetXaxis()->GetXmin();
        auto xmax = input_histo[bin_idx-1][lambda_idx][layer]->GetXaxis()->GetXmax();
        std::unique_ptr<TF1> fitfunc = std::make_unique<TF1>("fitfunc", "gaus", xmin, xmax);
        input_histo[bin_idx-1][lambda_idx][layer]->Fit("fitfunc", "Q");
        auto chi2 = fitfunc->GetChisquare();
        auto dof = fitfunc->GetNDF();
        
        if (dof)
        {
            double tmp_goodness = chi2/dof;
            if (goodness[layer] == 999)
            {
                goodness[layer] = tmp_goodness;
                best_lambda[layer] = lambda_start + lambda_step*lambda_idx;
                best_lambda_idx[layer] = lambda_idx;
            }
            else
            {
                if (abs(tmp_goodness-1) < abs(goodness[layer]-1))
                {
                    goodness[layer] = tmp_goodness;
                    best_lambda[layer] = lambda_start + lambda_step*lambda_idx;
                    best_lambda_idx[layer] = lambda_idx;
                }
            }
        }
    }
}