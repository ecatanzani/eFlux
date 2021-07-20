#include <vector>
#include <memory>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iostream>

#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TMath.h"
#include "TGraphErrors.h"

const double generation_vertex_radius = 1.381976597885342;
const int npx = 100000;
const int nGFitLoops = 10000;
const double chiSQLLimit = 0.000001;

std::vector<double> getEnergyBinning(const char* config_file);
std::vector<double> createLogBinning(const double eMin, const double eMax, const std::size_t n_bins);
double wtsydp(const double minene, const double maxene, const double index);
void tmpFit(std::shared_ptr<TF1> fitter_func, TGraphErrors &gr);
void tmpFit(TF1 &fitter_func, TH1D* histo);
void setStartingParameters(TF1 &oldFitter, TF1 &newFitter);

double logisticFunction(double *x, double *par);
double logisticFunction_st1(double *x, double *par);
double logisticFunction_st2(double *x, double *par);
double logisticFunction_st3(double *x, double *par);
double logisticFunction_st4(double *x, double *par);
double logisticFunction_st5(double *x, double *par);
double logisticFunction_st6(double *x, double *par);
double logisticFunction_st7(double *x, double *par);
double logisticFunction_st8(double *x, double *par);
double logisticFunction_st9(double *x, double *par);

void buildAcceptance(
    const char* acc_config_file,
    const char* full_acceptance_histos, 
    const char* output_file = "acceptance.root")
{
    auto energy_binning = getEnergyBinning(acc_config_file);

    TFile* input_file = TFile::Open(full_acceptance_histos, "READ");
    if (input_file->IsZombie())
    {
        std::cerr << "\n\nError opening input file: [" << full_acceptance_histos << "]\n\n";
        exit(100);
    }

    TH1D* h_gen = static_cast<TH1D*>(input_file->Get("h_gen"));
    TH1D* h_geometric = static_cast<TH1D*>(input_file->Get("h_geometric"));
    TH1D* h_geometric_trigger = static_cast<TH1D*>(input_file->Get("h_geometric_trigger"));
    TH1D* h_bgo_fiducial = static_cast<TH1D*>(input_file->Get("h_bgo_fiducial"));
    TH1D* h_all_cut = static_cast<TH1D*>(input_file->Get("h_all_cut"));

    h_gen->SetDirectory(0);
    h_geometric->SetDirectory(0);
    h_geometric_trigger->SetDirectory(0);
    h_bgo_fiducial->SetDirectory(0);
    h_all_cut->SetDirectory(0);

    input_file->Close();

    TH1D* h_geometric_factor = static_cast<TH1D*>(h_geometric->Clone("h_geometric_factor"));
    TH1D* h_acc_geometric_trigger = static_cast<TH1D*>(h_geometric_trigger->Clone("h_acc_geometric_trigger"));
    TH1D* h_acc_bgo_fiducial = static_cast<TH1D*>(h_bgo_fiducial->Clone("h_acc_bgo_fiducial"));
    TH1D* h_acc_all_cut = static_cast<TH1D*>(h_all_cut->Clone("h_acc_all_cut"));

    h_geometric_factor->Divide(h_gen);
    h_acc_geometric_trigger->Divide(h_gen);
    h_acc_bgo_fiducial->Divide(h_gen);
    h_acc_all_cut->Divide(h_gen);

    double genSurface = 4 * TMath::Pi() * pow(generation_vertex_radius, 2) / 2;
    double scaleFactor = TMath::Pi() * genSurface;

    h_geometric_factor->Scale(scaleFactor);
    h_acc_geometric_trigger->Scale(scaleFactor);
    h_acc_bgo_fiducial->Scale(scaleFactor);
    h_acc_all_cut->Scale(scaleFactor);

    h_geometric_factor->GetYaxis()->SetTitle("geometric factor [m^{2} sr]");
    h_acc_geometric_trigger->GetYaxis()->SetTitle("acceptance [m^{2} sr]");
    h_acc_bgo_fiducial->GetYaxis()->SetTitle("acceptance [m^{2} sr]");
    h_acc_all_cut->GetYaxis()->SetTitle("acceptance [m^{2} sr]");

    std::vector<double> geometric_factor (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_geometric_trigger (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_bgo_fiducial (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_all_cut (h_geometric_factor->GetNbinsX(), 0);

    std::vector<double> geometric_factor_err (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_geometric_trigger_err (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_bgo_fiducial_err (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> acc_all_cut_err (h_geometric_factor->GetNbinsX(), 0);

    std::vector<double> energy (h_geometric_factor->GetNbinsX(), 0);
    std::vector<double> energy_err (h_geometric_factor->GetNbinsX(), 0);

    for (unsigned int idx=0; idx<energy.size(); ++idx)
    {
        energy[idx] = wtsydp(energy_binning[idx], energy_binning[idx+1], -2.7);
        geometric_factor[idx] = h_geometric_factor->GetBinContent(idx+1);
        geometric_factor_err[idx] = h_geometric_factor->GetBinError(idx+1);
        acc_geometric_trigger[idx] = h_acc_geometric_trigger->GetBinContent(idx+1);
        acc_geometric_trigger_err[idx] = h_acc_geometric_trigger->GetBinError(idx+1);
        acc_bgo_fiducial[idx] = h_acc_bgo_fiducial->GetBinContent(idx+1);
        acc_bgo_fiducial_err[idx] = h_acc_bgo_fiducial->GetBinError(idx+1);
        acc_all_cut[idx] = h_acc_all_cut->GetBinContent(idx+1);
        acc_all_cut_err[idx] = h_acc_all_cut->GetBinError(idx+1);
    }
    
    TGraphErrors gr_geometric_factor((int)geometric_factor.size(), &energy[0], &geometric_factor[0], &energy_err[0], &geometric_factor_err[0]);
    TGraphErrors gr_acc_geometric_trigger((int)geometric_factor.size(), &energy[0], &acc_geometric_trigger[0], &energy_err[0], &acc_geometric_trigger_err[0]);
    TGraphErrors gr_acc_bgo_fiducial((int)geometric_factor.size(), &energy[0], &acc_bgo_fiducial[0], &energy_err[0], &acc_bgo_fiducial_err[0]);
    TGraphErrors gr_acc_all_cut((int)geometric_factor.size(), &energy[0], &acc_all_cut[0], &energy_err[0], &acc_all_cut_err[0]);

    gr_geometric_factor.SetName("gr_geometric_factor");
    gr_acc_geometric_trigger.SetName("gr_acc_geometric_trigger");
    gr_acc_bgo_fiducial.SetName("gr_acc_bgo_fiducial");
    gr_acc_all_cut.SetName("gr_acc_all_cut");

    gr_geometric_factor.SetTitle("Geometric Factor");
    gr_acc_geometric_trigger.SetTitle("Acceptance - geometric + trigger");
    gr_acc_bgo_fiducial.SetTitle("Acceptance - BGO fiducial");
    gr_acc_all_cut.SetTitle("Acceptance - all cuts");

    gr_geometric_factor.GetXaxis()->SetTitle("Energy [GeV]");
    gr_acc_geometric_trigger.GetXaxis()->SetTitle("Energy [GeV]");
    gr_acc_bgo_fiducial.GetXaxis()->SetTitle("Energy [GeV]");
    gr_acc_all_cut.GetXaxis()->SetTitle("Energy [GeV]");

    gr_geometric_factor.GetYaxis()->SetTitle("geometric factor [m^{2} sr]");
    gr_acc_geometric_trigger.GetYaxis()->SetTitle("acceptance [m^{2} sr]");
    gr_acc_bgo_fiducial.GetYaxis()->SetTitle("acceptance [m^{2} sr]");
    gr_acc_all_cut.GetYaxis()->SetTitle("acceptance [m^{2} sr]");

    double fit_func_min = energy[0];
    double fit_func_max = energy[(int)geometric_factor.size()-1];

    std::cout << "\nFitting acceptances in [" << fit_func_min << ":" << fit_func_max << "] GeV\n";
    
    std::shared_ptr<TF1> fitter_func_geometric_factor = std::make_shared<TF1>("fitter_func_geometric_factor", "pol1", 2, fit_func_max);
    gr_geometric_factor.Fit("fitter_func_geometric_factor", "IRQ");

    // First stage fit
    TF1 fitter_func_geometric_trigger_st0("fitter_func_geometric_trigger_st0", logisticFunction, fit_func_min, fit_func_max, 2);
    TF1 fitter_func_bgo_fiducial_st0("fitter_func_bgo_fiducial_st0", logisticFunction, fit_func_min, fit_func_max, 2);
    TF1 fitter_func_all_cuts_st0("fitter_func_all_cuts_st0", logisticFunction, fit_func_min, fit_func_max, 2);

    fitter_func_geometric_trigger_st0.SetLineColor(kOrange-3);
    fitter_func_bgo_fiducial_st0.SetLineColor(kOrange-3);
    fitter_func_all_cuts_st0.SetLineColor(kOrange-3);
    
    fitter_func_geometric_trigger_st0.SetNpx(npx);
    fitter_func_bgo_fiducial_st0.SetNpx(npx);
    fitter_func_all_cuts_st0.SetNpx(npx);
    
    fitter_func_geometric_trigger_st0.SetParameters(10, 10);
    fitter_func_bgo_fiducial_st0.SetParameters(1, 1);
    fitter_func_all_cuts_st0.SetParameters(1, 1);

    tmpFit(fitter_func_geometric_trigger_st0, h_acc_geometric_trigger);
    tmpFit(fitter_func_bgo_fiducial_st0, h_acc_bgo_fiducial);
    tmpFit(fitter_func_all_cuts_st0, h_acc_all_cut);

    // Second stage fit
    TF1 fitter_func_geometric_trigger_st1("fitter_func_geometric_trigger_st1", logisticFunction_st1, fit_func_min, fit_func_max, 4);
    TF1 fitter_func_bgo_fiducial_st1("fitter_func_bgo_fiducial_st1", logisticFunction_st1, fit_func_min, fit_func_max, 4);
    TF1 fitter_func_all_cuts_st1("fitter_func_all_cuts_st1", logisticFunction_st1, fit_func_min, fit_func_max, 4);

    setStartingParameters(fitter_func_geometric_trigger_st0, fitter_func_geometric_trigger_st1);
    setStartingParameters(fitter_func_bgo_fiducial_st0, fitter_func_bgo_fiducial_st1);
    setStartingParameters(fitter_func_all_cuts_st0, fitter_func_all_cuts_st1);

    fitter_func_geometric_trigger_st1.SetLineColor(kRed+2);
    fitter_func_bgo_fiducial_st1.SetLineColor(kRed+2);
    fitter_func_all_cuts_st1.SetLineColor(kRed+2);
    
    fitter_func_geometric_trigger_st1.SetNpx(npx);
    fitter_func_bgo_fiducial_st1.SetNpx(npx);
    fitter_func_all_cuts_st1.SetNpx(npx);

    tmpFit(fitter_func_geometric_trigger_st1, h_acc_geometric_trigger);
    tmpFit(fitter_func_bgo_fiducial_st1, h_acc_bgo_fiducial);
    tmpFit(fitter_func_all_cuts_st1, h_acc_all_cut);

    // Third stage fit
    TF1 fitter_func_geometric_trigger_st2("fitter_func_geometric_trigger_st2", logisticFunction_st2, fit_func_min, fit_func_max, 6);
    TF1 fitter_func_bgo_fiducial_st2("fitter_func_bgo_fiducial_st2", logisticFunction_st2, fit_func_min, fit_func_max, 6);
    TF1 fitter_func_all_cuts_st2("fitter_func_all_cuts_st2", logisticFunction_st2, fit_func_min, fit_func_max, 6);

    setStartingParameters(fitter_func_geometric_trigger_st1, fitter_func_geometric_trigger_st2);
    setStartingParameters(fitter_func_bgo_fiducial_st1, fitter_func_bgo_fiducial_st2);
    setStartingParameters(fitter_func_all_cuts_st1, fitter_func_all_cuts_st2);

    fitter_func_geometric_trigger_st2.SetLineColor(kCyan+3);
    fitter_func_bgo_fiducial_st2.SetLineColor(kCyan+3);
    fitter_func_all_cuts_st2.SetLineColor(kCyan+3);
    
    fitter_func_geometric_trigger_st2.SetNpx(npx);
    fitter_func_bgo_fiducial_st2.SetNpx(npx);
    fitter_func_all_cuts_st2.SetNpx(npx);

    tmpFit(fitter_func_geometric_trigger_st2, h_acc_geometric_trigger);
    tmpFit(fitter_func_bgo_fiducial_st2, h_acc_bgo_fiducial);
    tmpFit(fitter_func_all_cuts_st2, h_acc_all_cut);

    // Fourth stage fit
    TF1 fitter_func_geometric_trigger_st3("fitter_func_geometric_trigger_st3", logisticFunction_st3, fit_func_min, fit_func_max, 7);
    TF1 fitter_func_bgo_fiducial_st3("fitter_func_bgo_fiducial_st3", logisticFunction_st3, fit_func_min, fit_func_max, 7);
    TF1 fitter_func_all_cuts_st3("fitter_func_all_cuts_st3", logisticFunction_st3, fit_func_min, fit_func_max, 7);

    setStartingParameters(fitter_func_geometric_trigger_st2, fitter_func_geometric_trigger_st3);
    setStartingParameters(fitter_func_bgo_fiducial_st2, fitter_func_bgo_fiducial_st3);
    setStartingParameters(fitter_func_all_cuts_st2, fitter_func_all_cuts_st3);

    fitter_func_geometric_trigger_st3.SetLineColor(kGreen+1);
    fitter_func_bgo_fiducial_st3.SetLineColor(kGreen+1);
    fitter_func_all_cuts_st3.SetLineColor(kGreen+1);
    
    fitter_func_geometric_trigger_st3.SetNpx(npx);
    fitter_func_bgo_fiducial_st3.SetNpx(npx);
    fitter_func_all_cuts_st3.SetNpx(npx);

    tmpFit(fitter_func_geometric_trigger_st3, h_acc_geometric_trigger);
    tmpFit(fitter_func_bgo_fiducial_st3, h_acc_bgo_fiducial);
    tmpFit(fitter_func_all_cuts_st3, h_acc_all_cut);

    // Fifth stage fit
    TF1 fitter_func_geometric_trigger_st4("fitter_func_geometric_trigger_st4", logisticFunction_st4, fit_func_min, fit_func_max, 8);
    TF1 fitter_func_bgo_fiducial_st4("fitter_func_bgo_fiducial_st4", logisticFunction_st4, fit_func_min, fit_func_max, 8);
    TF1 fitter_func_all_cuts_st4("fitter_func_all_cuts_st4", logisticFunction_st4, fit_func_min, fit_func_max, 8);

    setStartingParameters(fitter_func_geometric_trigger_st3, fitter_func_geometric_trigger_st4);
    setStartingParameters(fitter_func_bgo_fiducial_st3, fitter_func_bgo_fiducial_st4);
    setStartingParameters(fitter_func_all_cuts_st3, fitter_func_all_cuts_st4);

    fitter_func_geometric_trigger_st4.SetLineColor(kMagenta+1);
    fitter_func_bgo_fiducial_st4.SetLineColor(kMagenta+1);
    fitter_func_all_cuts_st4.SetLineColor(kMagenta+1);
    
    fitter_func_geometric_trigger_st4.SetNpx(npx);
    fitter_func_bgo_fiducial_st4.SetNpx(npx);
    fitter_func_all_cuts_st4.SetNpx(npx);

    tmpFit(fitter_func_geometric_trigger_st4, h_acc_geometric_trigger);
    tmpFit(fitter_func_bgo_fiducial_st4, h_acc_bgo_fiducial);
    tmpFit(fitter_func_all_cuts_st4, h_acc_all_cut);

    // Sixth stage fit
    TF1 fitter_func_geometric_trigger_st5("fitter_func_geometric_trigger_st5", logisticFunction_st5, fit_func_min, fit_func_max, 9);
    TF1 fitter_func_bgo_fiducial_st5("fitter_func_bgo_fiducial_st5", logisticFunction_st5, fit_func_min, fit_func_max, 9);
    TF1 fitter_func_all_cuts_st5("fitter_func_all_cuts_st5", logisticFunction_st5, fit_func_min, fit_func_max, 9);

    setStartingParameters(fitter_func_geometric_trigger_st4, fitter_func_geometric_trigger_st5);
    setStartingParameters(fitter_func_bgo_fiducial_st4, fitter_func_bgo_fiducial_st5);
    setStartingParameters(fitter_func_all_cuts_st4, fitter_func_all_cuts_st5);

    fitter_func_geometric_trigger_st5.SetLineColor(kBlue-7);
    fitter_func_bgo_fiducial_st5.SetLineColor(kBlue-7);
    fitter_func_all_cuts_st5.SetLineColor(kBlue-7);
    
    fitter_func_geometric_trigger_st5.SetNpx(npx);
    fitter_func_bgo_fiducial_st5.SetNpx(npx);
    fitter_func_all_cuts_st5.SetNpx(npx);

    tmpFit(fitter_func_geometric_trigger_st5, h_acc_geometric_trigger);
    tmpFit(fitter_func_bgo_fiducial_st5, h_acc_bgo_fiducial);
    tmpFit(fitter_func_all_cuts_st5, h_acc_all_cut);

    // Seventh stage fit
    TF1 fitter_func_geometric_trigger_st6("fitter_func_geometric_trigger_st6", logisticFunction_st6, fit_func_min, fit_func_max, 10);
    TF1 fitter_func_bgo_fiducial_st6("fitter_func_bgo_fiducial_st6", logisticFunction_st6, fit_func_min, fit_func_max, 10);
    TF1 fitter_func_all_cuts_st6("fitter_func_all_cuts_st6", logisticFunction_st6, fit_func_min, fit_func_max, 10);

    setStartingParameters(fitter_func_geometric_trigger_st5, fitter_func_geometric_trigger_st6);
    setStartingParameters(fitter_func_bgo_fiducial_st5, fitter_func_bgo_fiducial_st6);
    setStartingParameters(fitter_func_all_cuts_st5, fitter_func_all_cuts_st6);

    fitter_func_geometric_trigger_st6.SetLineColor(kRed-7);
    fitter_func_bgo_fiducial_st6.SetLineColor(kRed-7);
    fitter_func_all_cuts_st6.SetLineColor(kRed-7);
    
    fitter_func_geometric_trigger_st6.SetNpx(npx);
    fitter_func_bgo_fiducial_st6.SetNpx(npx);
    fitter_func_all_cuts_st6.SetNpx(npx);

    tmpFit(fitter_func_geometric_trigger_st6, h_acc_geometric_trigger);
    tmpFit(fitter_func_bgo_fiducial_st6, h_acc_bgo_fiducial);
    tmpFit(fitter_func_all_cuts_st6, h_acc_all_cut);

    // Eighth stage fit
    TF1 fitter_func_geometric_trigger_st7("fitter_func_geometric_trigger_st7", logisticFunction_st7, fit_func_min, fit_func_max, 12);
    TF1 fitter_func_bgo_fiducial_st7("fitter_func_bgo_fiducial_st7", logisticFunction_st7, fit_func_min, fit_func_max, 12);
    TF1 fitter_func_all_cuts_st7("fitter_func_all_cuts_st7", logisticFunction_st7, fit_func_min, fit_func_max, 12);

    setStartingParameters(fitter_func_geometric_trigger_st6, fitter_func_geometric_trigger_st7);
    setStartingParameters(fitter_func_bgo_fiducial_st6, fitter_func_bgo_fiducial_st7);
    setStartingParameters(fitter_func_all_cuts_st6, fitter_func_all_cuts_st7);

    fitter_func_geometric_trigger_st7.SetLineColor(kPink-7);
    fitter_func_bgo_fiducial_st7.SetLineColor(kPink-7);
    fitter_func_all_cuts_st7.SetLineColor(kPink-7);
    
    fitter_func_geometric_trigger_st7.SetNpx(npx);
    fitter_func_bgo_fiducial_st7.SetNpx(npx);
    fitter_func_all_cuts_st7.SetNpx(npx);

    tmpFit(fitter_func_geometric_trigger_st7, h_acc_geometric_trigger);
    tmpFit(fitter_func_bgo_fiducial_st7, h_acc_bgo_fiducial);
    tmpFit(fitter_func_all_cuts_st7, h_acc_all_cut);

    // Ninth stage fit
    TF1 fitter_func_geometric_trigger_st8("fitter_func_geometric_trigger_st8", logisticFunction_st8, fit_func_min, fit_func_max, 14);
    TF1 fitter_func_bgo_fiducial_st8("fitter_func_bgo_fiducial_st8", logisticFunction_st8, fit_func_min, fit_func_max, 14);
    TF1 fitter_func_all_cuts_st8("fitter_func_all_cuts_st8", logisticFunction_st8, fit_func_min, fit_func_max, 14);

    setStartingParameters(fitter_func_geometric_trigger_st7, fitter_func_geometric_trigger_st8);
    setStartingParameters(fitter_func_bgo_fiducial_st7, fitter_func_bgo_fiducial_st8);
    setStartingParameters(fitter_func_all_cuts_st7, fitter_func_all_cuts_st8);

    fitter_func_geometric_trigger_st8.SetLineColor(kGreen+4);
    fitter_func_bgo_fiducial_st8.SetLineColor(kGreen+4);
    fitter_func_all_cuts_st8.SetLineColor(kGreen+4);
    
    fitter_func_geometric_trigger_st8.SetNpx(npx);
    fitter_func_bgo_fiducial_st8.SetNpx(npx);
    fitter_func_all_cuts_st8.SetNpx(npx);

    tmpFit(fitter_func_geometric_trigger_st8, h_acc_geometric_trigger);
    tmpFit(fitter_func_bgo_fiducial_st8, h_acc_bgo_fiducial);
    tmpFit(fitter_func_all_cuts_st8, h_acc_all_cut);

    // Tenth stage fit
    TF1 fitter_func_geometric_trigger_st9("fitter_func_geometric_trigger_st9", logisticFunction_st9, fit_func_min, fit_func_max, 16);
    TF1 fitter_func_bgo_fiducial_st9("fitter_func_bgo_fiducial_st9", logisticFunction_st9, fit_func_min, fit_func_max, 16);
    TF1 fitter_func_all_cuts_st9("fitter_func_all_cuts_st9", logisticFunction_st9, fit_func_min, fit_func_max, 16);

    setStartingParameters(fitter_func_geometric_trigger_st8, fitter_func_geometric_trigger_st9);
    setStartingParameters(fitter_func_bgo_fiducial_st8, fitter_func_bgo_fiducial_st9);
    setStartingParameters(fitter_func_all_cuts_st8, fitter_func_all_cuts_st9);

    fitter_func_geometric_trigger_st9.SetLineColor(kCyan-3);
    fitter_func_bgo_fiducial_st9.SetLineColor(kCyan-3);
    fitter_func_all_cuts_st9.SetLineColor(kCyan-3);
    
    fitter_func_geometric_trigger_st9.SetNpx(npx);
    fitter_func_bgo_fiducial_st9.SetNpx(npx);
    fitter_func_all_cuts_st9.SetNpx(npx);

    tmpFit(fitter_func_geometric_trigger_st9, h_acc_geometric_trigger);
    tmpFit(fitter_func_bgo_fiducial_st9, h_acc_bgo_fiducial);
    tmpFit(fitter_func_all_cuts_st9, h_acc_all_cut);

    TFile* outfile = TFile::Open(output_file, "RECREATE");
    if (outfile->IsZombie())
    {
        std::cerr << "\n\nError writing output acceptance file: [" << output_file << "]\n\n";
        exit(100);
    }

    outfile->mkdir("GeometricFactor");
    outfile->cd("GeometricFactor");

    h_geometric_factor->Write();
    gr_geometric_factor.Write();
    fitter_func_geometric_factor->Write();

    outfile->mkdir("GeometricFactorTrigger");
    outfile->cd("GeometricFactorTrigger");

    h_acc_geometric_trigger->Write();
    gr_acc_geometric_trigger.Write();
    fitter_func_geometric_trigger_st0.Write();
    fitter_func_geometric_trigger_st1.Write();
    fitter_func_geometric_trigger_st2.Write();
    fitter_func_geometric_trigger_st3.Write();
    fitter_func_geometric_trigger_st4.Write();
    fitter_func_geometric_trigger_st5.Write();
    fitter_func_geometric_trigger_st6.Write();
    fitter_func_geometric_trigger_st7.Write();
    fitter_func_geometric_trigger_st8.Write();
    fitter_func_geometric_trigger_st9.Write();

    outfile->mkdir("BGOFiducialVolume");
    outfile->cd("BGOFiducialVolume");

    h_acc_bgo_fiducial->Write();
    gr_acc_bgo_fiducial.Write();
    fitter_func_bgo_fiducial_st0.Write();
    fitter_func_bgo_fiducial_st1.Write();
    fitter_func_bgo_fiducial_st2.Write();
    fitter_func_bgo_fiducial_st3.Write();
    fitter_func_bgo_fiducial_st4.Write();
    fitter_func_bgo_fiducial_st5.Write();
    fitter_func_bgo_fiducial_st6.Write();
    fitter_func_bgo_fiducial_st7.Write();
    fitter_func_bgo_fiducial_st8.Write();
    fitter_func_bgo_fiducial_st9.Write();

    outfile->mkdir("AllCuts");
    outfile->cd("AllCuts");

    h_acc_all_cut->Write();
    gr_acc_all_cut.Write();
    fitter_func_all_cuts_st0.Write();
    fitter_func_all_cuts_st1.Write();
    fitter_func_all_cuts_st2.Write();
    fitter_func_all_cuts_st3.Write();
    fitter_func_all_cuts_st4.Write();
    fitter_func_all_cuts_st5.Write();
    fitter_func_all_cuts_st6.Write();
    fitter_func_all_cuts_st7.Write();
    fitter_func_all_cuts_st8.Write();
    fitter_func_all_cuts_st9.Write();
    
    outfile->Close();

    std::cout << "\n\nOutput file has been written: [" << output_file << "]\n\n";
}

std::vector<double> getEnergyBinning(const char* config_file)
{
    std::size_t n_bins = 0;
    double min_event_energy = -999;
    double max_event_energy = -999;
    std::vector<double> energy_binning;

    std::ifstream input_file(config_file);
    if (!input_file.is_open())
	{
		std::cerr << "\nInput acceptance config file not found [" << config_file << "]\n\n";
		exit(100);
	}
    std::string input_string(
		(std::istreambuf_iterator<char>(input_file)),
		(std::istreambuf_iterator<char>()));
	input_file.close();
    std::string tmp_str;
	std::istringstream input_stream(input_string);
	std::string::size_type sz;

    while (input_stream >> tmp_str)
	{
		if (!strcmp(tmp_str.c_str(), "n_energy_bins"))
			input_stream >> n_bins;
		if (!strcmp(tmp_str.c_str(), "min_event_energy"))
		{
			input_stream >> tmp_str;
			min_event_energy = stod(tmp_str, &sz);
		}
		if (!strcmp(tmp_str.c_str(), "max_event_energy"))
		{
			input_stream >> tmp_str;
			max_event_energy = stod(tmp_str, &sz);
		}
	}

    return createLogBinning(min_event_energy, max_event_energy, n_bins);
}   

std::vector<double> createLogBinning(
    const double eMin,
	const double eMax,
	const std::size_t n_bins)
{
    std::vector<double> binning(n_bins+1, 0);
	double log_interval = (log10(eMax) - log10(eMin)) / n_bins;
	for (unsigned int bIdx = 0; bIdx <= n_bins; ++bIdx)
		binning[bIdx] = pow(10, log10(eMin) + bIdx * log_interval);

	return binning;
}

double wtsydp(const double minene, const double maxene, const double index)
{
    auto dene = maxene - minene;
    if (index != -1)
        return pow(fabs((pow(maxene, index + 1) - pow(minene, index + 1)) / ((index + 1) * dene)), 1. / index);
    else
        return dene / log(maxene / minene);
}

void setStartingParameters(TF1 &oldFitter, TF1 &newFitter)
{
    for (int pIdx = 0; pIdx < oldFitter.GetNpar(); ++pIdx)
        newFitter.SetParameter(pIdx, oldFitter.GetParameter(pIdx));
}

void tmpFit(std::shared_ptr<TF1> fitter_func, TGraphErrors &gr)
{
    auto nPars = fitter_func->GetNpar();
    double chisq = 0;
    
    for (int fIdx = 0; fIdx < nGFitLoops; ++fIdx)
    {
        if (!fIdx)
            gr.Fit(fitter_func->GetName(), "IREX0");
        else
        {
            for(auto pIdx = 0; pIdx < nPars; ++pIdx)
                fitter_func->SetParameter(pIdx, fitter_func->GetParameter(pIdx));
            gr.Fit(fitter_func->GetName(), "IREX0");
            
            if (fabs(chisq - fitter_func->GetChisquare()) < chiSQLLimit)
                break;
        }
        chisq = fitter_func->GetChisquare();
        std::cout << "\nchi2: " << chisq;
    }
}

void tmpFit(TF1 &fitter_func, TH1D* histo)
{
    auto nPars = fitter_func.GetNpar();
    double chisq = 0;

    for (int fIdx = 0; fIdx < nGFitLoops; ++fIdx)
    {
        if (!fIdx)
            histo->Fit(fitter_func.GetName(), "WMIRN");
        else
        {
            for(auto pIdx = 0; pIdx < nPars; ++pIdx)
                fitter_func.SetParameter(pIdx, fitter_func.GetParameter(pIdx));
            histo->Fit(fitter_func.GetName(), "WMIRN");
            if (fabs(chisq - fitter_func.GetChisquare()) < chiSQLLimit)
                break;
        }
        chisq = fitter_func.GetChisquare();
    }
}

double logisticFunction(double *x, double *par)
{
    double xx = log10(x[0]+1);
    double logistic = par[0] / (1 + TMath::Exp(-par[1]*xx));
    return logistic;
}

double logisticFunction_st1(double *x, double *par)
{
    double xx = log10(x[0]+1);
    double logistic = par[0] / (1 + TMath::Exp(-par[1] * xx + par[2]*TMath::Power(xx,2) + par[3]));
    return logistic;
}

double logisticFunction_st2(double *x, double *par)
{
    double xx = log10(x[0]+1);
    double logistic = (par[0] + par[4] * xx + par[5] * TMath::Power(xx,2) ) / (1 + TMath::Exp(-par[1] * xx + par[2]*TMath::Power(xx,2) + par[3]));
    return logistic;
}

double logisticFunction_st3(double *x, double *par)
{
    double xx = log10(x[0]+1);
    double logistic = (par[0] + par[4] * xx + par[5] * TMath::Power(xx,2)) / (1 + TMath::Exp(-par[1] * xx + par[2]*TMath::Power(xx,2) + par[3] + par[6]*TMath::Power(xx,3) ));
    return logistic;
}

double logisticFunction_st4(double *x, double *par)
{
    double xx = log10(x[0]+1);
    double logistic = (par[0] + par[4] * xx + par[5] * TMath::Power(xx,2)) / (1 + TMath::Exp(-par[1] * xx + par[2]*TMath::Power(xx,2) + par[3] + par[6]*TMath::Power(xx,3) + par[7]*TMath::Power(xx,4) ));
    return logistic;
}

double logisticFunction_st5(double *x, double *par)
{
    double xx = log10(x[0]+1);
    double logistic = (par[0] + par[4] * xx + par[5] * TMath::Power(xx,2) + par[8] * TMath::Power(xx, 3)) / (1 + TMath::Exp(-par[1] * xx + par[2]*TMath::Power(xx,2) + par[3] + par[6]*TMath::Power(xx,3) + par[7]*TMath::Power(xx,4) ));
    return logistic;
}

double logisticFunction_st6(double *x, double *par)
{
    double xx = log10(x[0]+1);
    double logistic = (par[0] + par[4] * xx + par[5] * TMath::Power(xx,2) + par[8] * TMath::Power(xx, 3) + par[9]*TMath::Power(xx, 4)) / (1 + TMath::Exp(-par[1] * xx + par[2]*TMath::Power(xx,2) + par[3] + par[6]*TMath::Power(xx,3) + par[7]*TMath::Power(xx,4) ));
    return logistic;
}

double logisticFunction_st7(double *x, double *par)
{
    double xx = log10(x[0]+1);
    double logistic = (par[0] + par[4] * xx + par[5] * TMath::Power(xx,2) + par[8] * TMath::Power(xx, 3) + par[9]*TMath::Power(xx, 4) + par[10]*TMath::Power(xx, 5) + par[11]*TMath::Power(xx, 7) ) / (1 + TMath::Exp(-par[1] * xx + par[2]*TMath::Power(xx,2) + par[3] + par[6]*TMath::Power(xx,3) + par[7]*TMath::Power(xx,4) ));
    return logistic;
}

double logisticFunction_st8(double *x, double *par)
{
    double xx = log10(x[0]+1);
    double logistic = (par[0] + par[4] * xx + par[5] * TMath::Power(xx,2) + par[8] * TMath::Power(xx, 3) + par[9]*TMath::Power(xx, 4) + par[10]*TMath::Power(xx, 5) + par[11]*TMath::Power(xx, 7) + par[12]*TMath::Power(xx, 9) + par[13]*TMath::Power(xx, 11)) / (1 + TMath::Exp(-par[1] * xx + par[2]*TMath::Power(xx,2) + par[3] + par[6]*TMath::Power(xx,3) + par[7]*TMath::Power(xx,4) ));
    return logistic;
}

double logisticFunction_st9(double *x, double *par)
{
    double xx = log10(x[0]+1);
    double logistic = (par[0] + par[4] * xx + par[5] * TMath::Power(xx,2) + par[8] * TMath::Power(xx, 3) + par[9]*TMath::Power(xx, 4) + par[10]*TMath::Power(xx, 5) + par[11]*TMath::Power(xx, 7) + par[12]*TMath::Power(xx, 9) + par[13]*TMath::Power(xx, 11) + par[14]*TMath::Power(xx, 13) + par[15]*TMath::Power(xx, 16)) / (1 + TMath::Exp(-par[1] * xx + par[2]*TMath::Power(xx,2) + par[3] + par[6]*TMath::Power(xx,3) + par[7]*TMath::Power(xx,4) ));
    return logistic;
}