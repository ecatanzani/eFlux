#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <stdio.h>  
#include <stdlib.h>
#include <math.h>       /* isnan, sqrt */ 

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TArrayD.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TMath.h"
#include "RVersion.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCategory.h"
#include "TMVA/IMethod.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
#include "TMVA/DataLoader.h"
#endif
#endif

#include "inc/RICHBDTEle.h"
#include "inc/trdmap.h"

using namespace std;
using namespace TMVA;

bool kBatch=false;

//int bin_start=54;
//int bin_stop=57;
// int bin_start=40;
// int bin_stop=57;
// int bin_start=30;
// int bin_stop=42;
// int bin_start=73;//>73 means >700GeV
// int bin_stop=76;
//mduranti@ui02-ams:BDTRICHEle> less ../TreeBuilder/rootuples_B1130_pass7_eos.txt | head -1
//root://eosams.cern.ch///eos/ams/Data/AMS02/2018/ISS.B1130/pass7/1305853512.00000001.root
//mduranti@ui02-ams:BDTRICHEle> less ../TreeBuilder/rootuples_B1130_pass7_eos.txt | tail -1
//root://eosams.cern.ch///eos/ams/Data/AMS02/2018/ISS.B1130/pass7/1572344387.00000001.root 
//mduranti@ui02-ams:BDTRICHEle> ls -lh TrainingAndTestTrees/ | head -2
//totale 44G
//-rw-r--r-- 1 mduranti ams  18K 21 mag 11.25 TrainingAndTestTrees_13055_v-1.root
//mduranti@ui02-ams:BDTRICHEle> ls -lh TrainingAndTestTrees/ | tail -1
//-rw-r--r-- 1 mduranti ams  18K 21 mag 14.05 TrainingAndTestTrees_15724_v-1.root
int bin_start=13055;
int bin_stop=15724;

int SCversion=-1; //VERSIONE DELLA BDT

TString TrainingAndTestTreesPath = "TrainingAndTestTrees";

/*
const int nbozzo=62;
double bozzo_bins[nbozzo]={
  0.50, 0.72, 0.98, 1.29, 1.65, 2.05, 2.51, 3.01, 3.57, 4.19,
  4.87, 5.60, 6.40, 7.26, 8.20, 9.20, 10.28, 11.44, 12.68, 14.01,
  15.43, 16.95, 18.57, 20.30, 22.14, 24.10, 26.18, 28.40, 30.76, 33.26,
  35.93, 38.76, 41.76, 44.96, 48.35, 51.94, 55.76, 59.81, 64.11, 68.67,
  73.51, 78.64, 84.08, 89.86, 95.98, 102.47, 109.36, 116.66, 124.40, 132.60,
  141.31, 150.53, 159.58, 170.55, 184.25, 202.06, 226.59, 263.57, 350.00, 500.00,
  700.00, 1000.00
};
*/

/*
const int nbozzo=59; //era questo
Double_t bozzo_bins[nbozzo] = {
  0.50,   0.72,   0.98,   1.29,   1.65,   2.05,   2.51,    3.01,     3.57,  4.19,
  4.87,   5.60,   6.40,   7.26,   8.20,   9.20,  10.28,   11.44,    12.68,  14.01,
  15.43,  16.95,  18.57,  20.30,  22.14,  24.10,  26.18,   28.40,    30.76,  33.26,
  35.93,  38.76,  41.76,  44.96,  48.35,  51.94,  55.76,   59.81,    64.11,  68.67,
  73.51,  78.64,  84.08,  89.86,  95.98, 102.47, 113.95,  126.57,   140.47, 155.76,
  172.59, 191.10, 220.00, 260.00, 350.00, 500.00, 700.00, 1000.00, 10000.00
};
*/

/*
const int nbozzo=77; //BINNING DELLO SPLITTING 
Double_t bozzo_bins[nbozzo] = {
  0.50, 0.65, 0.82, 1.01, 1.22,
  1.46, 1.72, 2.00, 2.31, 2.65,
  3.00, 3.36, 3.73, 4.12, 4.54,
  5.00, 5.49, 6.00, 6.54, 7.10,
  7.69, 8.30, 8.95, 9.62, 10.32,
  11.04, 11.80, 12.59, 13.41, 14.25,
  15.14, 16.05, 17.00, 17.98, 18.99,
  20.04, 21.13, 22.25, 23.42, 24.62,
  25.90, 27.25, 28.68, 30.21, 31.82,
  33.53, 35.36, 37.31, 39.39, 41.61,
  44.00, 46.57, 49.33, 52.33, 55.58,
  59.13, 63.02, 67.30, 72.05, 77.37,
  83.36, 90.19, 98.08, 107.34, 118.42,
  132.11, 148.81, 169.86, 197.69, 237.16,
  290.0, 370.0, 500.0, 700.0, 1000.0,
  1500.00, 2200.00
};
*/

const int nbozzo=99999;
Double_t bozzo_bins[nbozzo];

int main(int argc, char *argv[]);
void SetBranchAddressAndActivate(TTree* t3, const char* bname, void* ptr);
void BuildTrees(bool kFast=false, int kDataOrMC=0);
void Train(bool kTestWithData=false);
void Plot();
TH2* SliceNormalizeX(TH2* h);
TH2* SliceNormalizeY(TH2* h);
TH1* CreateHistoFromGraph(TGraph* gr);

//---------------------------------------------------------

int main(int argc, char *argv[]){

  kBatch=true;
  
  TString argv0 = argv[0];

  if (argc>=3){
    bin_start=atoi(argv[1]);
    bin_stop=atoi(argv[2]);
  }

  printf("------------------------------------------\n");
  printf("Processing bins [%d, %d]\n", bin_start, bin_stop);
  printf("------------------------------------------\n");

  if (argv0.Contains("Train")) {
    if (argc==4) TrainingAndTestTreesPath=argv[3];

    if (argv0.Contains("TestWithData")) Train(true); 
    else Train();
  }
  else if (argv0.Contains("BuildTrees")) {
    if (argc==4) SCversion=atoi(argv[3]);

    if (SCversion<-1 || SCversion>1) { printf("SCversion:%d not valid\n",SCversion); return 7; }
    else if (SCversion==-1) printf("Not filling the BDT branch\n");
    
    if (argv0.Contains("Fast")) BuildTrees(true); 
    else BuildTrees();
  }
  else if (argv0.Contains("Plot")) {
    if (argc==4) TrainingAndTestTreesPath=argv[3];
    else TrainingAndTestTreesPath = "PlotTrees";

    SCversion=0;//the files, unfortunately, has this in the name so we need to know...
    if (argc==5) SCversion=atoi(argv[4]);
    if (SCversion<-1 || SCversion>1) { printf("SCversion:%d not valid\n", SCversion); return 7; }

    Plot();
  }
  else {
    printf("%s is the generic name of the executable. You have to call one of the various BuildTrees, Train, Plots, etc...\n", argv0.Data());
    return -1;
  }
  
  return 0;
}

//---------------------------------------------------------

void SetBranchAddressAndActivate(TTree* t3, const char* bname, void* ptr){

  t3->SetBranchStatus(bname, 1);
  t3->SetBranchAddress(bname, ptr);

  return;
}

//---------------------------------------------------------

//void BuildTrees
// opens Data and MC superskimmed files
// For every energy bin writes one output file (writes data+mc in the same file)
// Selects good events and saves them in separate trees (signal/background test/train)
// Saves only relevant variables
//     +If you already have BDT weights writes also BDT in the tree
void BuildTrees(bool kFast, int kDataOrMC) {
  //kDataOrMC:
  // - 0 everything from data
  
  TString opt = "RECREATE";
  //  if (kDataOrMC==1 || kDataOrMC==2 || kDataOrMC==3) opt = "UPDATE"; 

  TFile* output[nbozzo]; //divide gli output per bin energia
  int bozzo_bin_start=bin_start;
  int bozzo_bin_stop = bin_stop;
  printf("bozzo_bin_start=%d, bozzo_bin_stop=%d\n", bozzo_bin_start, bozzo_bin_stop);
  for (int ii=bozzo_bin_start; ii<=bozzo_bin_stop; ii++) {
    printf("init output[%d]\n", ii);
    output[ii] = new TFile(Form("TrainingAndTestTrees_%02d_v%d.root", ii, SCversion), opt.Data());
    printf("    --> %s\n", output[ii]->GetName());
  }

  //MC
  TTree* signal[nbozzo]={0};
  TTree* signal_test[nbozzo]={0};
  TTree* background[nbozzo]={0};
  TTree* background_test[nbozzo]={0};
  //Dati
  TTree* signal_test_dati[nbozzo]={0};
  TTree* background_test_dati[nbozzo]={0};

  //----------------------------------------------------

  //Declaration of leaves types
  UInt_t          Run;
  //   UInt_t          RunTag;
  UInt_t          Event;
  UInt_t          UTime;       //HeaderR::UTime 
  UInt_t          Time;        //JLV1 time (sec)
  Int_t           runAnalysisTag;

  Float_t         trdKlhrEP;
  Float_t         trdKlhrHP;
  Float_t         trdKlhrEH;
  // Float_t         trdKlhrEP_ene;
  // Float_t         trdKlhrHP_ene;
  // Float_t         trdKlhrEH_ene;
  // Float_t         trdKlhE_ene;
  // Float_t         trdKlhP_ene;
  // Float_t         trdKlhH_ene;
  // Float_t         trdKlhrEP_ene_trdstd;
  // Float_t         trdKlhrHP_ene_trdstd;
  // Float_t         trdKlhrEH_ene_trdstd;
  // Float_t         trdKlhE_ene_trdstd;
  // Float_t         trdKlhP_ene_trdstd;
  // Float_t         trdKlhH_ene_trdstd;
  Short_t         trdKHitsOnTrkTrack;
  Short_t         trdKHitsOnTrdTrack;
  
  Float_t         ecalBDTv5Es;
  Float_t         ecalBDTv6Es;
  Float_t         ecalBDTv7Es;
  Float_t         ecalCorrEneE_E;
  Float_t         ecalCorrEneE_P;
  Float_t         ecalCorrEneP_E;
  Float_t         ecalCorrEneP_P;
  
  Short_t         nTrTracks;
  Short_t         trkSpan;
  Bool_t          trkHasLay2;
  Float_t         trkDefRig;
  Float_t         trkDefTheta;                        //TrTrackR::GetTheta( default )
  Float_t         trkDefPhi;                          //TrTrackR::GetPhi( default )
  Float_t         trkInnQ;
    
  Int_t           trigPhysBPatt;

  Float_t         tofBeta;
  Float_t         tofQ;
  Float_t         tofUpperQ;
  Float_t         tofLowerQ;
  Float_t         tofQlay[4];

  int             nRichRing;
  float           richBeta;
  float           richBetaError;
  float           richProb;
  bool            richGoodBeta;
  bool            richIsGood;
  bool            richIsClean;
  float           richMass;
  bool            richIsNaF;
  float           richExpectedPhotoElectrons_true;
  float           richExpectedPhotoElectrons_false;
  float           richPhotoElectrons_true;
  float           richPhotoElectrons_false;
 
  Float_t         rtiLifetime;
  Float_t         rtiMaxCutoff[4][2];
  Float_t         rtiMaxIGRFCutoff[4][2];
  Int_t           passRTI;
  
  //----------------------------------------------------

  TChain *ch[1];//only on data for now
  Long64_t total_entries=0; 

  for (int dd=0; dd<1; dd++) {
    
    ch[dd] = new TChain("tree");
    
    int ntrees=0;
  
    //----------------------------------------------------
    //---Mac---
    //      TString pdata = "/Volumes/CaseSensitive2/ISS.B950.pass6.21Dec2015.Extended15Jan2017.20Jan2017.SuperskimmEne6.Hadded/";
    //---CNAF---
    TString pdata = "/storage/gpfs_ams/ams/groups/PGAllElectrons/ISS.B1130.pass7.01May2020.ECALOrRICH/";
    //----------------------------------------------------
    for (int ii=bin_start; ii<=bin_stop; ii++) {
      // ntrees += ch[dd]->Add(Form("%s/Positive_ForAnal_%02d*.root", pdata.Data(), ii));
      // ntrees += ch[dd]->Add(Form("%s/Negative_ForAnal_%02d*.root", pdata.Data(), ii));
      // ntrees += ch->Add(Form("%s/Positive_NotForAnal*%02d*.root", pdata.Data(), ii));
      // ntrees += ch->Add(Form("%s/Negative_NotForAnal*%02d*.root", pdata.Data(), ii));
      ntrees += ch[dd]->Add(Form("%s/%d*_skimmed.root", pdata.Data(), ii));
    }
  
    ch[dd]->GetListOfFiles()->Print();
  
    Long64_t nentries = ch[dd]->GetEntries();
    printf("TChain with %lld entries (%d trees)\n", nentries, ntrees);
    //  sleep(5);
  
    ch[dd]->SetBranchStatus("*",0);  // disable all branches
  
    // Set branch addresses.
    SetBranchAddressAndActivate(ch[dd], "Run",           &Run);
    SetBranchAddressAndActivate(ch[dd], "Event",         &Event);
    SetBranchAddressAndActivate(ch[dd], "Time",          &Time);
    SetBranchAddressAndActivate(ch[dd], "UTime",         &UTime);
    SetBranchAddressAndActivate(ch[dd], "runAnalysisTag", &runAnalysisTag);
  
    SetBranchAddressAndActivate(ch[dd], "trdKlhrEP", &trdKlhrEP);
    SetBranchAddressAndActivate(ch[dd], "trdKlhrHP", &trdKlhrHP);
    SetBranchAddressAndActivate(ch[dd], "trdKlhrEH", &trdKlhrEH);
    // SetBranchAddressAndActivate(ch[dd], "trdKlhrEP_ene", &trdKlhrEP_ene);
    // SetBranchAddressAndActivate(ch[dd], "trdKlhrHP_ene", &trdKlhrHP_ene);
    // SetBranchAddressAndActivate(ch[dd], "trdKlhrEH_ene", &trdKlhrEH_ene);
    // SetBranchAddressAndActivate(ch[dd], "trdKlhE_ene", &trdKlhE_ene);
    // SetBranchAddressAndActivate(ch[dd], "trdKlhP_ene", &trdKlhP_ene);
    // SetBranchAddressAndActivate(ch[dd], "trdKlhH_ene", &trdKlhH_ene);
    // SetBranchAddressAndActivate(ch[dd], "trdKlhrEP_ene_trdstd", &trdKlhrEP_ene_trdstd);
    // SetBranchAddressAndActivate(ch[dd], "trdKlhrHP_ene_trdstd", &trdKlhrHP_ene_trdstd);
    // SetBranchAddressAndActivate(ch[dd], "trdKlhrEH_ene_trdstd", &trdKlhrEH_ene_trdstd);
    // SetBranchAddressAndActivate(ch[dd], "trdKlhE_ene_trdstd", &trdKlhE_ene_trdstd);
    // SetBranchAddressAndActivate(ch[dd], "trdKlhP_ene_trdstd", &trdKlhP_ene_trdstd);
    // SetBranchAddressAndActivate(ch[dd], "trdKlhH_ene_trdstd", &trdKlhH_ene_trdstd);
    SetBranchAddressAndActivate(ch[dd], "trdKHitsOnTrkTrack", &trdKHitsOnTrkTrack);
    SetBranchAddressAndActivate(ch[dd], "trdKHitsOnTrdTrack", &trdKHitsOnTrdTrack);
    
    SetBranchAddressAndActivate(ch[dd], "ecalBDTv5Es", &ecalBDTv5Es);
    SetBranchAddressAndActivate(ch[dd], "ecalBDTv6Es", &ecalBDTv6Es);
    SetBranchAddressAndActivate(ch[dd], "ecalBDTv7Es", &ecalBDTv7Es);
    SetBranchAddressAndActivate(ch[dd], "ecalCorrEneE_E",&ecalCorrEneE_E);
    SetBranchAddressAndActivate(ch[dd], "ecalCorrEneE_P", &ecalCorrEneE_P);
    SetBranchAddressAndActivate(ch[dd], "ecalCorrEneP_E", &ecalCorrEneP_E);
    SetBranchAddressAndActivate(ch[dd], "ecalCorrEneP_P", &ecalCorrEneP_P);
    
    SetBranchAddressAndActivate(ch[dd], "nTrTracks", &nTrTracks);
    SetBranchAddressAndActivate(ch[dd], "trkSpan", &trkSpan);
    SetBranchAddressAndActivate(ch[dd], "trkHasLay2", &trkHasLay2);
    SetBranchAddressAndActivate(ch[dd], "trkDefRig", &trkDefRig);
    SetBranchAddressAndActivate(ch[dd], "trkDefTheta", &trkDefTheta);
    SetBranchAddressAndActivate(ch[dd], "trkDefPhi", &trkDefPhi);
    SetBranchAddressAndActivate(ch[dd], "trkInnQ", &trkInnQ);
    
    SetBranchAddressAndActivate(ch[dd], "trigPhysBPatt", &trigPhysBPatt);

    SetBranchAddressAndActivate(ch[dd], "tofBeta", &tofBeta);
    SetBranchAddressAndActivate(ch[dd], "tofQ", &tofQ);
    SetBranchAddressAndActivate(ch[dd], "tofUpperQ", &tofUpperQ);
    SetBranchAddressAndActivate(ch[dd], "tofLowerQ", &tofLowerQ);
    SetBranchAddressAndActivate(ch[dd], "tofQlay", tofQlay);

    SetBranchAddressAndActivate(ch[dd], "nRichRing",       &nRichRing);
    SetBranchAddressAndActivate(ch[dd], "richBeta",        &richBeta);
    SetBranchAddressAndActivate(ch[dd], "richBetaError",   &richBetaError);
    SetBranchAddressAndActivate(ch[dd], "richProb",        &richProb);
    SetBranchAddressAndActivate(ch[dd], "richGoodBeta",    &richGoodBeta);
    SetBranchAddressAndActivate(ch[dd], "richIsGood",      &richIsGood);
    SetBranchAddressAndActivate(ch[dd], "richIsClean",     &richIsClean);
    SetBranchAddressAndActivate(ch[dd], "richMass",        &richMass);
    SetBranchAddressAndActivate(ch[dd], "richIsNaF",       &richIsNaF);
    SetBranchAddressAndActivate(ch[dd], "richExpectedPhotoElectrons_true",        &richExpectedPhotoElectrons_true);
    SetBranchAddressAndActivate(ch[dd], "richExpectedPhotoElectrons_false",       &richExpectedPhotoElectrons_false);
    SetBranchAddressAndActivate(ch[dd], "richPhotoElectrons_true",                &richPhotoElectrons_true);
    SetBranchAddressAndActivate(ch[dd], "richPhotoElectrons_false",               &richPhotoElectrons_false);
    
    if (dd==0) { //Data only
      SetBranchAddressAndActivate(ch[dd],  "rtiLifetime", &rtiLifetime);
      SetBranchAddressAndActivate(ch[dd],  "rtiMaxCutoff", rtiMaxCutoff);
      SetBranchAddressAndActivate(ch[dd],  "rtiMaxIGRFCutoff", rtiMaxIGRFCutoff);
      SetBranchAddressAndActivate(ch[dd],  "passRTI",&passRTI);
    }
    
    //----------------------------------------------------------------
    printf("All branches set\n");
    //----------------------------------------------------------------

    for (int ii=bozzo_bin_start; ii<=bozzo_bin_stop; ii++) {
      printf("cd-ing in %s...   ii=%d  dd=%d\n", output[ii]->GetName(), ii, dd);
      output[ii]->cd();
      signal[ii]               = ch[dd]->CloneTree(0);
      signal_test[ii]          = ch[dd]->CloneTree(0);
      signal[ii]->SetName("TreeS");//Signal
      signal_test[ii]->SetName("TreeST");//SignalTest
      background[ii]           = ch[dd]->CloneTree(0);
      background_test[ii]      = ch[dd]->CloneTree(0);
      background[ii]->SetName("TreeB");//Background
      background_test[ii]->SetName("TreeBT");//BackgroundTest
    }
  
    //----------------------------------------------------------------
    printf("All output trees created\n");
    //----------------------------------------------------------------
  
    double BDT;
    double BDTG;
  
    for (int ii=bozzo_bin_start; ii<=bozzo_bin_stop; ii++) {
      TTree* p[6] = {signal[ii], signal_test[ii], background[ii], background_test[ii], signal_test_dati[ii], background_test_dati[ii]};
    
      int start=-99;
      int stop=-99;
    
      if (dd==0) {
	start=0;
	stop=3;
      }
    
      for (int jj=start; jj<=stop; jj++) {
	p[jj]->Branch("BDT",  &BDT,  "BDT/D");
	p[jj]->Branch("BDTG", &BDTG, "BDTG/D");
      }
    }
  
    //---------------------------------------------------------------
  
    int perc=1;
  
    Long64_t nbytes = 0;
  
    TTree *_signal=NULL;
    TTree *_signal_test=NULL;
    TTree *_background=NULL;
    TTree *_background_test=NULL;
    TTree *_signal_test_dati=NULL;
    TTree *_background_test_dati=NULL;
  
    bool kQUICK=false;
    //    double kPrescale = 0.3;
    double kPrescale = 0.04;
    double kAdditionalProtonPrescale = 0.01;
    if (SCversion>=0) { 
      kPrescale = 1.0;
      kAdditionalProtonPrescale = 1.0;
    }

    printf("Prescale = %f, AdditionalProtonPrescale = %f\n", kPrescale, kAdditionalProtonPrescale);
    
    Long64_t ee=0;
    for (ee=0; ee<nentries; ee++) {

      if (gRandom->Uniform(0,1)>kPrescale) continue;
      
      total_entries++; 
    
      nbytes += ch[dd]->GetEntry(ee);
      double pperc=(100.0*ee)/nentries;
      if (pperc>perc) {
	printf("Processed %lld out of %lld: %d%%\n", ee+1, nentries, (int)perc);
	perc++;
      }
      //      printf("Processing %lld...\n", ee);
      //      printf("Event=%u, Run=%u: passed std-cuts\n", Event, Run);
      //-------------------------------------------------------------------------------------------------------------------
      if (kQUICK && ee>100 ) break; 
    
      _signal = NULL;
      _signal_test = NULL;
      _background = NULL;
      _background_test = NULL;
      _signal_test_dati = NULL;
      _background_test_dati = NULL;
    
      BDT = -999;
      BDTG = -999;
    
      double richBeta_capped = std::min((double)richBeta, 1.0);
      richMass = fabs(trkDefRig/richBeta_capped)*sqrt(1-richBeta_capped*richBeta_capped);

      //-------------------------------------------------------------------------------------------------------------------
    
      bool cFisici = (trigPhysBPatt & 0x3E)!=0;
      bool cTrdHits = trdKHitsOnTrkTrack>8;
      bool cTrdDefined = trdKlhrEP>0;
      //      bool cTrdDefined = trdKlhE_ene>0;
      bool cCharge = trkInnQ<1.5;
      bool cTrdeHe = trdKlhrEH<0.8;
      //      bool cTrdeHe = trdKlhrEH_ene<0.8;
      bool cSingleTrack = nTrTracks<2;
      bool cGood = runAnalysisTag>=0 && passRTI>0 && rtiLifetime>0;
      //      bool cCutOff = ecalCorrEneE_E>1.25*TMath::Max(fabs(rtiMaxCutoff[0][0]), fabs(rtiMaxCutoff[0][1]));
      //      bool cCutOff = fabs(trkDefRig)>1.25*TMath::Max(fabs(rtiMaxCutoff[0][0]), fabs(rtiMaxCutoff[0][1]));
      bool cCutOff = fabs(trkDefRig)>TMath::Max(fabs(rtiMaxIGRFCutoff[0][0]), fabs(rtiMaxIGRFCutoff[0][1]));
    
      bool cTrue = 1;
    
      bool stdcuts =
	cFisici &&
	cTrdHits &&
	cTrdDefined &&
	cCharge &&
	cTrdeHe &&
	cSingleTrack && 
	cGood &&
	//	cCutOff &&
	cTrue;
    
      //-------------------------------------------------------------------------------------------------------------------
    
      bool cNegative = trkDefRig<0.0;
      bool cPositive = trkDefRig>0.0;

      //-------------------------------------------------------------------------------------------------------------------

      bool cRICHCut = nRichRing>0 && richBeta>0;
      bool cECALCut = ecalCorrEneE_E>0.5;

      //-------------------------------------------------------------------------------------------------------------------

      // attenzione questo cut sull'energia non si puo' togliere perche' per rari eventi potrei avere degli arrotondamenti per alcuni eventi
      //      bool cAddCut = ecalCorrEneE_E>=bozzo_bins[bozzo_bin_start] && ecalCorrEneE_E<bozzo_bins[bozzo_bin_stop];
      //      bool cAddCut = ecalEnergyD>=bozzo_bins[bozzo_bin_start] && ecalEnergyD<bozzo_bins[bozzo_bin_stop]; 
      bool cAddCut = cRICHCut && fabs(trkDefRig)<25.0; //this is required for the trees to train the RICH BDT. To produce the trees for the plots the one below is used
      if (SCversion>=0) {//this is used to produce the trees, with all the pre-selected events, for the plots
	cAddCut = fabs(trkDefRig)<25.0;
      }
      
      //      float trdclass = -(log10(trdKlhE_ene)+2);
    
      //--------------------------------------------------------------------------------
      // bool cElectron = (trdKlhrEP_ene<0.5 && trkDefRig<0);
      // bool cProton =  (trdKlhrEP_ene>0.9 && trkDefRig>0); // ceci: trd killer cut
      //--------------------------------------------------------------------------------
      // bool cElectron = (trdKlhrEP_ene<0.6 && trkDefRig<0);
      // bool cProton =  (trdKlhrEP_ene>0.75 && trkDefRig>0);
      //--------------------------------------------------------------------------------
      // bool cElectron = (trdclass>1.0 && cNegative);
      // bool cProton = (trdclass<1.0 && cPositive);
      //--------------------------------------------------------------------------------
      bool cElectron = (trdKlhrEP<0.5 && cNegative);
      bool cProton =  (trdKlhrEP>0.8 && cPositive);
      //--------------------------------------------------------------------------------        

      if (SCversion>=0) {
	cElectron = cNegative;
	cProton =  cPositive;
      }
      
      // printf("Energy = %f\n", ecalCorrEneE_E);
    
      if (stdcuts) {
	
	if (cAddCut) {//this is used to produce the trees, with all the pre-selected events, for the plots
	  //printf("Event=%u, Run=%u: passed std-cuts\n", Event, Run);
	
	  //---------------------------------------------------------------
	
	  //for (int ii=bin_start; ii<=bin_stop; ii++) { // valid only if file splitting bin and bozzobins are equal
	  TString sRun = Form("%d", Run);
	  for (int ii=bozzo_bin_start; ii<=bozzo_bin_stop; ii++) { 
	    //	    if (ecalCorrEneE_E>=bozzo_bins[ii] && ecalCorrEneE_E<bozzo_bins[ii+1]) {
	    //	    if (ecalEnergyD>=bozzo_bins[ii] && ecalEnergyD<bozzo_bins[ii+1]) {
	    TString sBozzoBinStart = Form("%d", ii);
	    if (sRun.BeginsWith(sBozzoBinStart)) {
	      //assegna il puntatore del bin giusto
	      _signal = signal[ii];
	      _signal_test = signal_test[ii];
	      _background = background[ii];
	      _background_test = background_test[ii];
	      _signal_test_dati = signal_test_dati[ii];
	      _background_test_dati = background_test_dati[ii];
	      // printf(" energy %f    bin %d  singnal_ %p \n",ecalCorrEneE_E,ii,signal_ );
	      // printf("%f --> %d\n", ecalCorrEneE_E, ii);
	      break;
	    }
	  }
	
	  //---------------------------------------------------------------
	
	  if (!(cElectron || cProton)) continue;
	
	  //---------------------------------------------------------------
	
	  if (SCversion>=0 && cRICHCut && (tofBeta==tofBeta)) { //sometimes tofBeta is NaN

	    struct data_for_BDT dfb;

	    dfb.tofBeta = tofBeta;
	    dfb.tofQ = tofQ;
	    dfb.tofUpperQ = tofUpperQ;
	    dfb.tofLowerQ = tofLowerQ;
	    for (int cc=0; cc<4; cc++) {
	      dfb.tofQlay[cc] = tofQlay[cc];
	    }
	    dfb.nRichRing = nRichRing;
	    dfb.richGoodBeta = richGoodBeta;
	    dfb.richIsGood = richIsGood;
	    dfb.richIsClean = richIsClean;
	    dfb.richIsNaF = richIsNaF;
	    dfb.richBeta = richBeta;
	    dfb.richBetaError = richBetaError;
	    dfb.richProb = richProb;
	    dfb.richMass = richMass;
	    dfb.richExpectedPhotoElectrons_true = richExpectedPhotoElectrons_true;
	    dfb.richExpectedPhotoElectrons_false = richExpectedPhotoElectrons_false;
	    dfb.richPhotoElectrons_true = richPhotoElectrons_true;
	    dfb.richPhotoElectrons_false = richPhotoElectrons_false; 

	    BDT = GetBDT(dfb, "BDT", SCversion);
	    BDTG = GetBDT(dfb, "BDTG", SCversion);
	  }
	
	  //---------------------------------------------------------------
	
	  if (dd==0){
	    if (cElectron) {
	      // printf("Event=%u, Run=%u: Filling S...\n", Event, Run);
	      if (gRandom->Uniform()<0.5) {//50% of events
		_signal->Fill();
	      }
	      else {
		_signal_test->Fill();
	      }
	    }
	    else if (cProton) {
	      if (gRandom->Uniform(0,1)>kAdditionalProtonPrescale) continue;
	      // printf("Event=%u, Run=%u: Filling B...\n", Event, Run);
	      if (gRandom->Uniform()<0.5) {//50% of events
		_background->Fill();
	      }
	      else {
		_background_test->Fill();
	      }
	    }
	  }
	}	 
      }
    
    }//loop on entries
    printf("All entries looped...\n");
  
  }//loop on chains
  printf("All input chains processed...\n");
  
  for (int ii=bozzo_bin_start; ii<=bozzo_bin_stop; ii++) {

    printf("Bin=%d) Processed=%lld, Accepted(S,ST,B,BT)=%lld,%lld,%lld,%lld = %f,%f,%f,%f\n", ii, total_entries,
	   signal[ii]->GetEntries(), signal_test[ii]->GetEntries(), background[ii]->GetEntries(), background_test[ii]->GetEntries(),
	   signal[ii]->GetEntries()/((double)total_entries), signal_test[ii]->GetEntries()/((double)total_entries), background[ii]->GetEntries()/((double)total_entries), background_test[ii]->GetEntries()/((double)total_entries)
	   );
    output[ii]->Write(output[ii]->GetName(), TObject::kOverwrite);
    output[ii]->Close();
  }
  
  return;
}

void Train(bool kTestWithData) {
  // The explicit loading of the shared libTMVA is done in TMVAlogon.C, defined in .rootrc
  // if you use your private .rootrc, or run from a different directory, please copy the
  // corresponding lines from .rootrc
  
  // methods to be processed can be given as an argument; use format:
  //
  // mylinux~> root -l TMVAClassification.C
  //
    
  //---------------------------------------------------------------
  // This loads the library
  TMVA::Tools::Instance();
  
  // to get access to the GUI and all tmva macros
  // TString tmva_dir(TString(gRootDir) + "/tmva");
  // gROOT->SetMacroPath(tmva_dir + "/test/:" + gROOT->GetMacroPath() );
  gROOT->ProcessLine(".L ./TMVAGui.C");
  
  // Default MVA methods to be trained + tested
  std::map<std::string,int> Use;
  
  // --- Cut optimisation
  Use["Cuts"]            = 0;
  Use["CutsD"]           = 0;
  Use["CutsPCA"]         = 0;
  Use["CutsGA"]          = 0;
  Use["CutsSA"]          = 0;
  // 
  // --- 1-dimensional likelihood ("naive Bayes estimator")
  Use["Likelihood"]      = 0;
  Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
  Use["LikelihoodPCA"]   = 0; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
  Use["LikelihoodKDE"]   = 0;
  Use["LikelihoodMIX"]   = 0;
  //
  // --- Mutidimensional likelihood and Nearest-Neighbour methods
  Use["PDERS"]           = 0;
  Use["PDERSD"]          = 0;
  Use["PDERSPCA"]        = 0;
  Use["PDEFoam"]         = 0;
  Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
  Use["KNN"]             = 0; // k-nearest neighbour method
  //
  // --- Linear Discriminant Analysis
  Use["LD"]              = 0; // Linear Discriminant identical to Fisher
  Use["Fisher"]          = 0;
  Use["FisherG"]         = 0;
  Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
  Use["HMatrix"]         = 0;
  //
  // --- Function Discriminant analysis
  Use["FDA_GA"]          = 0; // minimisation of user-defined function using Genetics Algorithm
  Use["FDA_SA"]          = 0;
  Use["FDA_MC"]          = 0;
  Use["FDA_MT"]          = 0;
  Use["FDA_GAMT"]        = 0;
  Use["FDA_MCMT"]        = 0;
  //
  // --- Neural Networks (all are feed-forward Multilayer Perceptrons)
  Use["MLP"]             = 0; // Recommended ANN
  Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
  Use["MLPBNN"]          = 0; // Recommended ANN with BFGS training method and bayesian regulator
  Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
  Use["TMlpANN"]         = 0; // ROOT's own ANN
  //
  // --- Support Vector Machine 
  Use["SVM"]             = 0;
  // 
  // --- Boosted Decision Trees
  Use["BDT"]             = 1; // uses Adaptive Boost
  Use["BDTG"]            = 1; // uses Gradient Boost
  Use["BDTB"]            = 0; // uses Bagging
  Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
  Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting 
  // 
  // --- Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
  Use["RuleFit"]         = 0;
  // ---------------------------------------------------------------

  std::cout << std::endl;
  std::cout << "==> Start TMVA Classification" << std::endl;

  // --------------------------------------------------------------------------------------------------

  // --- Here the preparation phase begins

  // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
  TString outfileName;
  if (kTestWithData) {
    outfileName = "TMVA_data.root";
  }
  else {
    outfileName = "TMVA.root";
  }
  TFile* outputFile = TFile::Open(outfileName, "RECREATE");
  
  // Create the factory object. Later you can choose the methods
  // whose performance you'd like to investigate. The factory is 
  // the only TMVA object you have to interact with
  //
  // The first argument is the base of the name of all the
  // weightfiles in the directory weight/
  //
  // The second argument is the output file for the training results
  // All TMVA output can be suppressed by removing the "!" (not) in
  // front of the "Silent" argument in the option string
  TString base_name = "TMVAClassification";
  TMVA::Factory* factory = new TMVA::Factory(base_name.Data(), outputFile,
					     Form(
						  // "!V:!Silent:%s:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
						  // "!V:!Silent:%s:DrawProgressBar:Transformations=I;D;P;D:AnalysisType=Classification");
						  // "!V:!Silent:%sColor:DrawProgressBar:Transformations=I;D:AnalysisType=Classification" );
						  "V:!Silent:%s:DrawProgressBar:Transformations=I:AnalysisType=Classification"
						  , kBatch?"!Color":"Color"));
  
  // If you wish to modify default settings
  // (please check "src/Config.h" to see all available global options)
  //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
  //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";

  // Define the input variables that shall be used for the MVA training
  // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
  // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
  TMVA::DataLoader* loader = new TMVA::DataLoader("dataset");
#define _LOADER loader,
#else
  TMVA::Factory* loader = factory;
#define _LOADER
#endif
  
  loader->AddVariable("tofBeta", "tofBeta", "units", 'F');
  loader->AddVariable("tofQ", "tofQ", "units", 'F');
  loader->AddVariable("tofUpperQ", "tofUpperQ", "units", 'F');
  loader->AddVariable("tofLowerQ", "tofLowerQ", "units", 'F');
  for (int ii=0; ii<4; ii++) {
    loader->AddVariable(Form("tofQlay[%d]", ii), Form("tofQlay_%d", ii), "units", 'F');
  }
  loader->AddVariable("nRichRing", "nRichRing", "units", 'I');
  loader->AddVariable("richGoodBeta", "richGoodBeta", "units", 'I');
  loader->AddVariable("richIsGood", "richIsGood", "units", 'I');
  loader->AddVariable("richIsClean", "richIsClean", "units", 'I');
  loader->AddVariable("richIsNaF", "richIsNaF", "units", 'I');
  
  loader->AddVariable("richBeta", "richBeta", "units", 'F');
  loader->AddVariable("richBetaError", "richBetaError", "units", 'F');
  loader->AddVariable("richProb", "richProb", "units", 'F');
  loader->AddVariable("richMass", "richMass", "units", 'F');
  loader->AddVariable("richExpectedPhotoElectrons_true", "richExpectedPhotoElectrons_true", "units", 'F');
  loader->AddVariable("richExpectedPhotoElectrons_false", "richExpectedPhotoElectrons_false", "units", 'F');
  loader->AddVariable("richPhotoElectrons_true", "richPhotoElectrons_true", "units", 'F');
  loader->AddVariable("richPhotoElectrons_false", "richPhotoElectrons_false", "units", 'F');
    
  // You can add so-called "Spectator variables", which are not used in the MVA training,
  // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
  // input variables, the response values of all trained MVAs, and the spectator variables
  loader->AddSpectator("ecalBDTv5Es", "ecalBDTv5Es", "units", 'F');
  loader->AddSpectator("ecalBDTv6Es", "ecalBDTv6Es", "units", 'F');
  loader->AddSpectator("ecalBDTv7Es", "ecalBDTv7Es", "units", 'F');

  loader->AddSpectator("trkDefRig", "trkDefRig", "units", 'F');
  //  loader->AddSpectator("trkInnQ", "trkInnQ", "units", 'F');

  loader->AddSpectator("trdKlhrEP", "trdKlhrEP", "units", 'F');
  loader->AddSpectator("trdKlhrHP", "trdKlhrHP", "units", 'F');
  loader->AddSpectator("trdKlhrEH", "trdKlhrEH", "units", 'F');
  // loader->AddSpectator("trdKlhrEP_ene", "trdKlhrEP_ene", "units", 'F');
  // loader->AddSpectator("trdKlhrHP_ene", "trdKlhrHP_ene", "units", 'F');
  // loader->AddSpectator("trdKlhrEH_ene", "trdKlhrEH_ene", "units", 'F');
  // loader->AddSpectator("trdKlhE_ene", "trdKlhE_ene", "units", 'F');
  // loader->AddSpectator("trdKlhP_ene", "trdKlhP_ene", "units", 'F');
  // loader->AddSpectator("trdKlhH_ene", "trdKlhH_ene", "units", 'F');
  // loader->AddSpectator("trdKlhrEP_ene_trdstd", "trdKlhrEP_ene_trdstd", "units", 'F');
  // loader->AddSpectator("trdKlhrHP_ene_trdstd", "trdKlhrHP_ene_trdstd", "units", 'F');
  // loader->AddSpectator("trdKlhrEH_ene_trdstd", "trdKlhrEH_ene_trdstd", "units", 'F');
  // loader->AddSpectator("trdKlhE_ene_trdstd", "trdKlhE_ene_trdstd", "units", 'F');
  // loader->AddSpectator("trdKlhP_ene_trdstd", "trdKlhP_ene_trdstd", "units", 'F');
  // loader->AddSpectator("trdKlhH_ene_trdstd", "trdKlhH_ene_trdstd", "units", 'F');

  printf("---------------------------------------------------\n");

  //---------------------------------------------------

  // --- Register the training and test trees
  TChain* TreeS   = new TChain("TreeS");  
  TChain* TreeST  = new TChain("TreeST");
  TChain* TreeB   = new TChain("TreeB");
  TChain* TreeBT  = new TChain("TreeBT");
  TChain* TreeSTD = new TChain("TreeSTD");
  TChain* TreeBTD = new TChain("TreeBTD");
  
  // Read training and test data
  printf("*** Reading input files from %s/\n", TrainingAndTestTreesPath.Data());
  for (int ii=bin_start; ii<=bin_stop; ii++) {
    TFile* input = new TFile(Form("%s/TrainingAndTestTrees_%02d_v%d.root", TrainingAndTestTreesPath.Data(), ii, SCversion));
    std::cout << "--- TMVAClassification       : Using input file: " << input->GetName() << std::endl;   
    if (!input) {
      cout << "ERROR: could not open data file" << endl;
    }
    else {
      TreeS->Add(input->GetName());
      TreeST->Add(input->GetName());
      TreeB->Add(input->GetName());
      TreeBT->Add(input->GetName());
      TreeSTD->Add(input->GetName());
      TreeBTD->Add(input->GetName());
      input->Close();
      delete input;
    }
  }

  TTree* signal = TreeS;

  TTree* signal_test = NULL;
  if (kTestWithData) signal_test = TreeSTD;
  else signal_test = TreeST;
  
  TTree* background = TreeB;
  
  TTree* background_test = NULL;
  if (kTestWithData) background_test = TreeBTD;
  else background_test = TreeBT;  
  
  if (signal->GetEntries()<1) {
    cout << "signal entries < 1" << endl;
    exit(1);
  }
  if (signal_test->GetEntries()<1) {
    cout << "signal_test entries < 1" << endl;
    exit(1);
  }
  if (background->GetEntries()<1) {
    cout << "background entries < 1" << endl;
    exit(1);
  }
  if (background_test->GetEntries()<1) {
    cout << "background_test entries < 1" << endl;
    exit(1);
  }

  Int_t nTrainSig = signal->GetEntries();
  Int_t nTestSig = signal_test->GetEntries();
  Int_t nTrainBkg = background->GetEntries();
  Int_t nTestBkg = background_test->GetEntries();
  
  // global event weights per tree (see below for setting event-wise weights)
  Double_t signalWeight     = 1.0;
  Double_t backgroundWeight = 1.0;

  // --- Register the training and test trees
  loader->AddSignalTree(signal, signalWeight, TMVA::Types::kTraining);
  loader->AddSignalTree(signal_test, signalWeight, TMVA::Types::kTesting);
  loader->AddBackgroundTree(background, backgroundWeight , TMVA::Types::kTraining);
  loader->AddBackgroundTree(background_test, backgroundWeight , TMVA::Types::kTesting);
      
  // You can add an arbitrary number of signal or background trees

  // Set individual event weights (the variables must exist in the original TTree)
  //    for signal    : loader->SetSignalWeightExpression    ("weight1*weight2");
  //    for background: loader->SetBackgroundWeightExpression("weight1*weight2");
  // loader->SetSignalWeightExpression("weight_maura*norm");
  // loader->SetBackgroundWeightExpression("weight_maura*norm");

  // --- end of tree registration 

  // Apply additional cuts on the signal and background samples (can be different)
  TCut mycuts = "";
  TCut mycutb = "";

  //  TCut energyrange = Form("ecalCorrEneE_E>=%f && ecalCorrEneE_E<%f", bozzo_bins[bin_start], bozzo_bins[bin_stop]);
  //  TCut energyrange = Form("ecalEnergyD>=%f && ecalEnergyD<%f", bozzo_bins[bin_start], bozzo_bins[bin_stop]);
  TCut energyrange = "";
  mycuts+=energyrange;
  mycutb+=energyrange;

  /*
  if (kTestWithData) {
    TCut lowenergy = Form("ecalCorrEneE_E < %f", 102.47);
    mycuts+=lowenergy;
    mycutb+=lowenergy;
  }
  */

  TCut removenans = "!(TMath::IsNaN(tofQ) || !(TMath::Finite(tofQ)))";
  mycuts+=removenans;
  mycutb+=removenans;

  TCut removenans2 = "!(TMath::IsNaN(tofQlay[1]) || !(TMath::Finite(tofQlay[1])))";
  mycuts+=removenans2;
  mycutb+=removenans2;
  
  mycuts.Print();
  mycutb.Print();

  // Tell the loader how to use the training and testing events
  //
  // If no numbers of events are given, half of the events in the tree are used 
  // for training, and the other half for testing:
  //    loader->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );
  // To also specify the number of testing events, use:
  loader->PrepareTrainingAndTestTree(mycuts, mycutb,
				     //"nTrain_Signal=0:nTrain_Background=0:nTest_Signal=0:nTest_Background=0:SplitMode=Random:NormMode=NumEvents:!V");
  				     //"nTrain_Signal=0:nTrain_Background=0:nTest_Signal=0:nTest_Background=0:SplitMode=Random:NormMode=None:!V");
				     //"nTrain_Signal=1000000:nTrain_Background=1000000:nTest_Signal=0:nTest_Background=0:SplitMode=Random:NormMode=NumEvents:!V");
  				     "nTrain_Signal=1000000:nTrain_Background=1000000:nTest_Signal=0:nTest_Background=0:SplitMode=Random:NormMode=None:!V");
  
  // For an example of the category classifier usage, see: TMVAClassificationCategory
  TMVA::MethodCategory* mcategory = 0;

  // ---- Book MVA methods
  //
  // Please lookup the various method configuration options in the corresponding cxx files, eg:
  // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
  // It is possible to preset ranges in the option string in which the cut optimisation should be done:
  // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

  // Cut optimisation
  if (Use["Cuts"])
    factory->BookMethod( _LOADER TMVA::Types::kCuts, "Cuts",
			 "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );

  if (Use["CutsD"])
    factory->BookMethod( _LOADER TMVA::Types::kCuts, "CutsD",
			 "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate" );

  if (Use["CutsPCA"])
    factory->BookMethod( _LOADER TMVA::Types::kCuts, "CutsPCA",
			 "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA" );

  if (Use["CutsGA"])
    factory->BookMethod( _LOADER TMVA::Types::kCuts, "CutsGA",
			 "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" );

  if (Use["CutsSA"])
    factory->BookMethod( _LOADER TMVA::Types::kCuts, "CutsSA",
			 "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

  // Likelihood ("naive Bayes estimator")
  if (Use["Likelihood"])
    factory->BookMethod( _LOADER TMVA::Types::kLikelihood, "Likelihood",
			 "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" );

  // Decorrelated likelihood
  if (Use["LikelihoodD"])
    factory->BookMethod( _LOADER TMVA::Types::kLikelihood, "LikelihoodD",
			 "!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" );

  // PCA-transformed likelihood
  if (Use["LikelihoodPCA"])
    factory->BookMethod( _LOADER TMVA::Types::kLikelihood, "LikelihoodPCA",
			 "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA" ); 

  // Use a kernel density estimator to approximate the PDFs
  if (Use["LikelihoodKDE"])
    factory->BookMethod( _LOADER TMVA::Types::kLikelihood, "LikelihoodKDE",
			 "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50" ); 

  // Use a variable-dependent mix of splines and kernel density estimator
  if (Use["LikelihoodMIX"])
    factory->BookMethod( _LOADER TMVA::Types::kLikelihood, "LikelihoodMIX",
			 "!H:!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50" ); 

  // Test the multi-dimensional probability density estimator
  // here are the options strings for the MinMax and RMS methods, respectively:
  //      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );
  //      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );
  if (Use["PDERS"])
    factory->BookMethod( _LOADER TMVA::Types::kPDERS, "PDERS",
			 "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" );

  if (Use["PDERSD"])
    factory->BookMethod( _LOADER TMVA::Types::kPDERS, "PDERSD",
			 "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate" );

  if (Use["PDERSPCA"])
    factory->BookMethod( _LOADER TMVA::Types::kPDERS, "PDERSPCA",
			 "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA" );

  // Multi-dimensional likelihood estimator using self-adapting phase-space binning
  if (Use["PDEFoam"])
    factory->BookMethod( _LOADER TMVA::Types::kPDEFoam, "PDEFoam",
			 "!H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" );

  if (Use["PDEFoamBoost"])
    factory->BookMethod( _LOADER TMVA::Types::kPDEFoam, "PDEFoamBoost",
			 "!H:!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T" );

  // K-Nearest Neighbour classifier (KNN)
  if (Use["KNN"])
    factory->BookMethod( _LOADER TMVA::Types::kKNN, "KNN",
			 "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );

  // H-Matrix (chi2-squared) method
  if (Use["HMatrix"])
    factory->BookMethod( _LOADER TMVA::Types::kHMatrix, "HMatrix", "!H:!V:VarTransform=None" );

  // Linear discriminant (same as Fisher discriminant)
  if (Use["LD"])
    factory->BookMethod( _LOADER TMVA::Types::kLD, "LD", "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

  // Fisher discriminant (same as LD)
  if (Use["Fisher"])
    factory->BookMethod( _LOADER TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

  // Fisher with Gauss-transformed input variables
  if (Use["FisherG"])
    factory->BookMethod( _LOADER TMVA::Types::kFisher, "FisherG", "H:!V:VarTransform=Gauss" );

  // Composite classifier: ensemble (tree) of boosted Fisher classifiers
  if (Use["BoostedFisher"])
    factory->BookMethod( _LOADER TMVA::Types::kFisher, "BoostedFisher", 
			 "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2:!Boost_DetailedMonitoring" );

  // Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
  if (Use["FDA_MC"])
    factory->BookMethod( _LOADER TMVA::Types::kFDA, "FDA_MC",
			 "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1" );

  if (Use["FDA_GA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
    factory->BookMethod( _LOADER TMVA::Types::kFDA, "FDA_GA",
			 "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1" );

  if (Use["FDA_SA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
    factory->BookMethod( _LOADER TMVA::Types::kFDA, "FDA_SA",
			 "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

  if (Use["FDA_MT"])
    factory->BookMethod( _LOADER TMVA::Types::kFDA, "FDA_MT",
			 "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );

  if (Use["FDA_GAMT"])
    factory->BookMethod( _LOADER TMVA::Types::kFDA, "FDA_GAMT",
			 "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );

  if (Use["FDA_MCMT"])
    factory->BookMethod( _LOADER TMVA::Types::kFDA, "FDA_MCMT",
			 "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20" );

  // TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
  if (Use["MLP"])
    factory->BookMethod( _LOADER TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );

  if (Use["MLPBFGS"])
    factory->BookMethod( _LOADER TMVA::Types::kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );

  if (Use["MLPBNN"])
    factory->BookMethod( _LOADER TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" ); // BFGS training with bayesian regulators

  // CF(Clermont-Ferrand)ANN
  if (Use["CFMlpANN"])
    factory->BookMethod( _LOADER TMVA::Types::kCFMlpANN, "CFMlpANN", "!H:!V:NCycles=2000:HiddenLayers=N+1,N"  ); // n_cycles:#nodes:#nodes:...  

  // Tmlp(Root)ANN
  if (Use["TMlpANN"])
    factory->BookMethod( _LOADER TMVA::Types::kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  ); // n_cycles:#nodes:#nodes:...

  // Support Vector Machine
  if (Use["SVM"])
    factory->BookMethod( _LOADER TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );

  // Boosted Decision Trees
  if (Use["BDTG"]) {// Gradient Boost
    factory->BookMethod( _LOADER TMVA::Types::kBDT, "BDTG",
			 // "!H:!V:NTrees=1000:BoostType=Grad:nCuts=200:NegWeightTreatment=NoNegWeightsInTraining");
    			 //"!H:!V:NTrees=100:BoostType=Grad:nCuts=200:NegWeightTreatment=NoNegWeightsInTraining:MaxDepth=10:PruneMethod=ExpectedError:PruneStrength=50:DoBoostMonitor=True:DoPreselection=True");
    "!H:!V:NTrees=100:BoostType=Grad:nCuts=200:NegWeightTreatment=NoNegWeightsInTraining:MaxDepth=10:PruneMethod=ExpectedError:PruneStrength=50:DoBoostMonitor=True");

    /*    
    //Categories in energies
    TMVA::MethodBase* bdtg_category = factory->BookMethod( _LOADER TMVA::Types::kCategory, "BDTG_Category", "" );
    mcategory = dynamic_cast<TMVA::MethodCategory*>(bdtg_category);
    
    for (int ii=bin_start; ii<=bin_stop; ii+=3){
      //     printf("ii = %d\n", ii);
      double first=bozzo_bins[ii];
      double last;
      if ((ii+3)<=bin_stop) last=bozzo_bins[ii+3];
      else last=bozzo_bins[bin_stop];
      printf("Creating category for %f<E<%f\n", first, last);
      TCut cut = Form("ecalEnergyE>=%f && ecalEnergyE<%f", first, last);
      cut.Print();
      mcategory->AddMethod( cut,
			    "F1SL:F2SL:F3SL:", TMVA::Types::kBDT, "Category_BDTG",
			    "!H:!V:NTrees=100:BoostType=Grad:UseBaggedGrad:nCuts=200" );
    }
    */

  }
  
  if (Use["BDT"])  // Adaptive Boost
    factory->BookMethod( _LOADER TMVA::Types::kBDT, "BDT",
			 // "!H:!V:NTrees=1000:BoostType=AdaBoost:nCuts=200" );
			 //			 "!H:!V:NTrees=100:BoostType=AdaBoost:nCuts=200:MaxDepth=10:PruneMethod=ExpectedError:PruneStrength=50:DoBoostMonitor=True:DoPreselection=True");
			 "!H:!V:NTrees=100:BoostType=AdaBoost:nCuts=200:MaxDepth=10:PruneMethod=ExpectedError:PruneStrength=50:DoBoostMonitor=True");
  
  if (Use["BDTB"]) // Bagging
    factory->BookMethod( _LOADER TMVA::Types::kBDT, "BDTB",
			 "!H:!V:NTrees=1000:BoostType=Bagging:nCuts=200" );
  
  if (Use["BDTD"]) {// Decorrelation + Adaptive Boost
    factory->BookMethod( _LOADER TMVA::Types::kBDT, "BDTD",
			 "!H:!V:NTrees=1000:BoostType=AdaBoost:nCuts=200:VarTransform=Decorrelate" );
    
    /*
    //Categories in energies
    TMVA::MethodBase* bdtd_category = factory->BookMethod( _LOADER TMVA::Types::kCategory, "BDTD_Category", "" );
    mcategory = dynamic_cast<TMVA::MethodCategory*>(bdtd_category);
    
    for (int ii=bin_start; ii<=bin_stop; ii+=3){
      //     printf("ii = %d\n", ii);
      double first=bozzo_bins[ii];
      double last;
      if ((ii+3)<=bin_stop) last=bozzo_bins[ii+3];
      else last=bozzo_bins[bin_stop];
      printf("Creating category for %f<E<%f\n", first, last);
      TCut cut = Form("ecalEnergyE>=%f && ecalEnergyE<%f", first, last);
      cut.Print();
      mcategory->AddMethod( cut,
			    "F1SL:F2SL:F3SL:", TMVA::Types::kBDT, "Category_BDTD",
			    "!H:!V:NTrees=100:BoostType=AdaBoost:nCuts=200:VarTransform=Decorrelate" );
    }
    */

  }

  if (Use["BDTF"])  // Allow Using Fisher discriminant in node splitting for (strong) linearly correlated variables
    factory->BookMethod( _LOADER TMVA::Types::kBDT, "BDTMitFisher",
			 "!H:!V:NTrees=1000:UseFisherCuts:BoostType=AdaBoost:nCuts=200" );

  // RuleFit -- TMVA implementation of Friedman's method
  if (Use["RuleFit"])
    factory->BookMethod( _LOADER TMVA::Types::kRuleFit, "RuleFit",
			 "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );
  
  // --------------------------------------------------------------------------------------------------
  
  // ---- Now you can optimize the setting (configuration) of the MVAs using the set of training events
  
  // loader->OptimizeAllMethods("SigEffAt001","Scan");
  // loader->OptimizeAllMethods("ROCIntegral","GA");
  
  // --------------------------------------------------------------------------------------------------
  
  // ---- Now you can tell the loader to train, test, and evaluate the MVAs

  // Train MVAs using the set of training events
  factory->TrainAllMethods();

  // ---- Evaluate all MVAs using the set of test events
  factory->TestAllMethods();

  // ----- Evaluate and compare performance of all configured MVAs
  factory->EvaluateAllMethods();

  // --------------------------------------------------------------

  // Save the output
  outputFile->Close();

  std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;

  delete loader;

#if defined(__CINT__) && not defined(__MAKECINT__)
  if (!gROOT->IsBatch()) TMVAGui( outfileName );
#endif

  return;
}

void Plot() {
  
  TFile* ftrdmap = new TFile("trdmap.root");
  trdmap* _trdmap = NULL;
  if (ftrdmap) {
    _trdmap = (trdmap*)(ftrdmap->Get("trdmap"));
  }
  if (_trdmap) {
    printf("*********************************************************\n");
    printf("trdmap with %d unique indeces and %d indeces loaded...\n", _trdmap->getmaxuniqueindex(), _trdmap->getmaxindex());
    printf("*********************************************************\n");
  }
  
  // Create a ROOT output file for histograms
  TString outfileName;
  outfileName = Form("Plots_%d_%d.root", bin_start, bin_stop);
  
  TFile* outputFile = TFile::Open(outfileName, "RECREATE");

  //----------------------------------------------------

  double start= 1305800000;
  double stop = 1572500000;
  int days = ((int)((stop-start)/(24.0*60.0*60.0)));

  TH1F* hEventsAll = new TH1F("hEventsAll", "All Events", days, start, stop);
  TH1F* hEventsRICHorECAL = new TH1F("hEventsRICHorECAL", "RICH or ECAL Events", days, start, stop);
  TH1F* hEventsECAL = new TH1F("hEventsECAL", "ECAL Events", days, start, stop);
  TH1F* hEventsRICH = new TH1F("hEventsRICH", "RICH Events", days, start, stop);
  TH1F* hEventsElectrons = new TH1F("hEventsElectrons", "Electrons Events", days, start, stop);
  TH1F* hEventsProtons = new TH1F("hEventsProtons", "Protons Events", days, start, stop);
  TH1F* hEventsElectronsECAL = new TH1F("hEventsElectronsECAL", "Electrons ECAL Events", days, start, stop);
  TH1F* hEventsElectronsRICH = new TH1F("hEventsElectronsRICH", "Electrons RICH Events", days, start, stop);
  TH1F* hEventsProtonsECAL = new TH1F("hEventsProtonsECAL", "Protons ECAL Events", days, start, stop);
  TH1F* hEventsProtonsRICH = new TH1F("hEventsProtonsRICH", "Protons RICH Events", days, start, stop);

  TH1F* hTRDAll = new TH1F("hTRDAll", "All TRD", 100, 0.0, 2.0);
  TH1F* hTRDECAL = new TH1F("hTRDECAL", "ECAL TRD", 100, 0.0, 2.0);
  TH1F* hTRDRICH = new TH1F("hTRDRICH", "RICH TRD", 100, 0.0, 2.0);
  TH1F* hTRDElectrons = new TH1F("hTRDElectrons", "Electrons TRD", 100, 0.0, 2.0);
  TH1F* hTRDProtons = new TH1F("hTRDProtons", "Protons TRD", 100, 0.0, 2.0);
  TH1F* hTRDProtonsCC = new TH1F("hTRDProtonsCC", "Protons (CC) TRD", 100, 0.0, 2.0);
  TH1F* hTRDProtonsCCHS = new TH1F("hTRDProtonsCCHS", "Protons (CC - High Statistics) TRD", 100, 0.0, 2.0);
  TH1F* hTRDElectronsECAL = new TH1F("hTRDElectronsECAL", "Electrons ECAL TRD", 100, 0.0, 2.0);
  TH1F* hTRDElectronsRICH = new TH1F("hTRDElectronsRICH", "Electrons RICH TRD", 100, 0.0, 2.0);
  TH1F* hTRDProtonsECAL = new TH1F("hTRDProtonsECAL", "Protons ECAL TRD", 100, 0.0, 2.0);
  TH1F* hTRDProtonsRICH = new TH1F("hTRDProtonsRICH", "Protons RICH TRD", 100, 0.0, 2.0);
  TH1F* hTRDProtonsCCECAL = new TH1F("hTRDProtonsCCECAL", "Protons (CC) ECAL TRD", 100, 0.0, 2.0);
  TH1F* hTRDProtonsCCRICH = new TH1F("hTRDProtonsCCRICH", "Protons (CC) RICH TRD", 100, 0.0, 2.0);
  TH1F* hTRDProtonsCCHSECAL = new TH1F("hTRDProtonsCCHSECAL", "Protons (CC - High Statistics) ECAL TRD", 100, 0.0, 2.0);
  TH1F* hTRDProtonsCCHSRICH = new TH1F("hTRDProtonsCCHSRICH", "Protons (CC - High Statistics) RICH TRD", 100, 0.0, 2.0);

  TH2F* hTRDvsRElectrons = new TH2F("hTRDvsRElectrons", "Electrons TRD vs log10(R)", 100, -1, 2, 100, 0.0, 2.0);
  TH2F* hTRDvsRElectronsECAL = new TH2F("hTRDvsRElectronsECAL", "Electrons ECAL TRD vs log10(R)", 100, -1, 2, 100, 0.0, 2.0);
  TH2F* hTRDvsRElectronsRICH = new TH2F("hTRDvsRElectronsRICH", "Electrons RICH TRD vs log10(R)", 100, -1, 2, 100, 0.0, 2.0);

  TH2F* hTRDvsRProtons = new TH2F("hTRDvsRProtons", "Protons TRD vs log10(R)", 100, -1, 2, 100, 0.0, 2.0);
  TH2F* hTRDvsRProtonsECAL = new TH2F("hTRDvsRProtonsECAL", "Protons ECAL TRD vs log10(R)", 100, -1, 2, 100, 0.0, 2.0);
  TH2F* hTRDvsRProtonsRICH = new TH2F("hTRDvsRProtonsRICH", "Protons RICH TRD vs log10(R)", 100, -1, 2, 100, 0.0, 2.0);

  TH2F* hTRDvsRProtonsCC = new TH2F("hTRDvsRProtonsCC", "Protons (CC) TRD vs log10(R)", 100, -1, 2, 100, 0.0, 2.0);
  TH2F* hTRDvsRProtonsCCECAL = new TH2F("hTRDvsRProtonsCCECAL", "Protons (CC) ECAL TRD vs log10(R)", 100, -1, 2, 100, 0.0, 2.0);
  TH2F* hTRDvsRProtonsCCRICH = new TH2F("hTRDvsRProtonsCCRICH", "Protons (CC) RICH TRD vs log10(R)", 100, -1, 2, 100, 0.0, 2.0);
  TH2F* hTRDvsRProtonsCCHS = new TH2F("hTRDvsRProtonsCCHS", "Protons (CC - High Statistics) TRD vs log10(R)", 100, -1, 2, 100, 0.0, 2.0);
  TH2F* hTRDvsRProtonsCCHSECAL = new TH2F("hTRDvsRProtonsCCHSECAL", "Protons (CC - High Statistics) ECAL TRD vs log10(R)", 100, -1, 2, 100, 0.0, 2.0);
  TH2F* hTRDvsRProtonsCCHSRICH = new TH2F("hTRDvsRProtonsCCHSRICH", "Protons (CC - High Statistics) RICH TRD vs log10(R)", 100, -1, 2, 100, 0.0, 2.0);
  
  TH2F* hBDTvsRNegative = new TH2F("hBDTvsRNegative", "BDT (R<0) vs log10(R)", 100, -1, 2, 100, -1.0, 1.0);
  TH2F* hBDTvsRPositive = new TH2F("hBDTvsRPositive", "BDT (R>0) vs log10(R)", 100, -1, 2, 100, -1.0, 1.0);
  TH2F* hBDTvsRElectrons = new TH2F("hBDTvsRElectrons", "BDT (R<0) vs log10(R)", 100, -1, 2, 100, -1.0, 1.0);
  TH2F* hBDTvsRProtonsCC = new TH2F("hBDTvsRProtonsCC", "BDT (R<0) vs log10(R)", 100, -1, 2, 100, -1.0, 1.0);
  TH2F* hBDTvsRProtons = new TH2F("hBDTvsRProtons", "BDT (R>0) vs log10(R)", 100, -1, 2, 100, -1.0, 1.0);

  TH2F* hBDTGvsRNegative = new TH2F("hBDTGvsRNegative", "BDTG (R<0) vs log10(R)", 100, -1, 2, 100, -1.0, 1.0);
  TH2F* hBDTGvsRPositive = new TH2F("hBDTGvsRPositive", "BDTG (R>0) vs log10(R)", 100, -1, 2, 100, -1.0, 1.0);
  TH2F* hBDTGvsRElectrons = new TH2F("hBDTGvsRElectrons", "BDTG (R<0) vs log10(R)", 100, -1, 2, 100, -1.0, 1.0);
  TH2F* hBDTGvsRProtonsCC = new TH2F("hBDTGvsRProtonsCC", "BDTG (R<0) vs log10(R)", 100, -1, 2, 100, -1.0, 1.0);
  TH2F* hBDTGvsRProtons = new TH2F("hBDTGvsRProtons", "BDTG (R>0) vs log10(R)", 100, -1, 2, 100, -1.0, 1.0);

  
  TH1F* hTRDElectrons0005 = new TH1F("hTRDElectrons0005", "TRD Electrons #theta=(0#circ-5#circ)", 100, 0.0, 2.0);
  TH1F* hTRDElectrons0510 = new TH1F("hTRDElectrons0510", "TRD Electrons #theta=(5#circ-10#circ)", 100, 0.0, 2.0);
  TH1F* hTRDElectrons1015 = new TH1F("hTRDElectrons1015", "TRD Electrons #theta=(10#circ-15#circ)", 100, 0.0, 2.0);
  TH1F* hTRDElectrons1520 = new TH1F("hTRDElectrons1520", "TRD Electrons #theta=(15#circ-20#circ)", 100, 0.0, 2.0);
  TH1F* hTRDElectrons2025 = new TH1F("hTRDElectrons2025", "TRD Electrons #theta=(20#circ-25#circ)", 100, 0.0, 2.0);
  TH1F* hTRDElectrons2530 = new TH1F("hTRDElectrons2530", "TRD Electrons #theta=(25#circ-30#circ)", 100, 0.0, 2.0);
  TH1F* hTRDElectrons3035 = new TH1F("hTRDElectrons3035", "TRD Electrons #theta=(30#circ-35#circ)", 100, 0.0, 2.0);

  TH1F* hTRDProtons0005 = new TH1F("hTRDProtons0005", "TRD Protons #theta=(0#circ-5#circ)", 100, 0.0, 2.0);
  TH1F* hTRDProtons0510 = new TH1F("hTRDProtons0510", "TRD Protons #theta=(5#circ-10#circ)", 100, 0.0, 2.0);
  TH1F* hTRDProtons1015 = new TH1F("hTRDProtons1015", "TRD Protons #theta=(10#circ-15#circ)", 100, 0.0, 2.0);
  TH1F* hTRDProtons1520 = new TH1F("hTRDProtons1520", "TRD Protons #theta=(15#circ-20#circ)", 100, 0.0, 2.0);
  TH1F* hTRDProtons2025 = new TH1F("hTRDProtons2025", "TRD Protons #theta=(20#circ-25#circ)", 100, 0.0, 2.0);
  TH1F* hTRDProtons2530 = new TH1F("hTRDProtons2530", "TRD Protons #theta=(25#circ-30#circ)", 100, 0.0, 2.0);
  TH1F* hTRDProtons3035 = new TH1F("hTRDProtons3035", "TRD Protons #theta=(30#circ-35#circ)", 100, 0.0, 2.0);

  TH2F* hTRDProtonsVsTime = new TH2F("hTRDProtonsVsTime", "Protons TRD Vs Time", days, start, stop, 100, 0.0, 2.0);
  TH2F* hTRDNormalizedProtonsVsTime = new TH2F("hTRDNormalizedProtonsVsTime", "Protons TRD Vs Time", days, start, stop, 100, 0.0, 2.0);
  TH2F* hTRDProtonsCCVsTime = new TH2F("hTRDProtonsCCVsTime", "Protons (CC) TRD Vs Time", days, start, stop, 100, 0.0, 2.0);
  TH2F* hTRDProtonsCCHSVsTime = new TH2F("hTRDProtonsCCHSVsTime", "Protons (CC - High Statistics) TRD Vs Time", days, start, stop, 100, 0.0, 2.0);
  TH2F* hTRDNormalizedProtonsCCVsTime = new TH2F("hTRDNormalizedProtonsCCVsTime", "Protons (CC) TRD Vs Time", days, start, stop, 100, 0.0, 2.0);
  TH2F* hTRDNormalizedProtonsCCHSVsTime = new TH2F("hTRDNormalizedProtonsCCHSVsTime", "Protons (CC - High Statistics) TRD Vs Time", days, start, stop, 100, 0.0, 2.0);
  TH2F* hTRDElectronsVsTime = new TH2F("hTRDElectronsVsTime", "Electrons TRD Vs Time", days, start, stop, 100, 0.0, 2.0);
  TH2F* hTRDNormalizedElectronsVsTime = new TH2F("hTRDNormalizedElectronsVsTime", "Electrons TRD Vs Time", days, start, stop, 100, 0.0, 2.0);

  //----------------------------------------------------
    
  //Declaration of leaves types
  UInt_t          Run;
  //   UInt_t          RunTag;
  UInt_t          Event;
  UInt_t          UTime;       //HeaderR::UTime 
  UInt_t          Time;        //JLV1 time (sec)
  Int_t           runAnalysisTag;

  Float_t         trdKlhrEP;
  Float_t         trdKlhrHP;
  Float_t         trdKlhrEH;
  // Float_t         trdKlhrEP_ene;
  // Float_t         trdKlhrHP_ene;
  // Float_t         trdKlhrEH_ene;
  // Float_t         trdKlhE_ene;
  // Float_t         trdKlhP_ene;
  // Float_t         trdKlhH_ene;
  // Float_t         trdKlhrEP_ene_trdstd;
  // Float_t         trdKlhrHP_ene_trdstd;
  // Float_t         trdKlhrEH_ene_trdstd;
  // Float_t         trdKlhE_ene_trdstd;
  // Float_t         trdKlhP_ene_trdstd;
  // Float_t         trdKlhH_ene_trdstd;
  Short_t         trdKHitsOnTrkTrack;
  Short_t         trdKHitsOnTrdTrack;
  
  Float_t         ecalBDTv5Es;
  Float_t         ecalBDTv6Es;
  Float_t         ecalBDTv7Es;
  Float_t         ecalCorrEneE_E;
  Float_t         ecalCorrEneE_P;
  Float_t         ecalCorrEneP_E;
  Float_t         ecalCorrEneP_P;
  
  Short_t         nTrTracks;
  Short_t         trkSpan;
  Bool_t          trkHasLay2;
  Float_t         trkDefRig;
  Float_t         trkDefTheta;                        //TrTrackR::GetTheta( default )
  Float_t         trkDefPhi;                          //TrTrackR::GetPhi( default )
  Float_t         trkInnQ;
    
  Int_t           trigPhysBPatt;

  Float_t         tofBeta;
  Float_t         tofQ;
  Float_t         tofUpperQ;
  Float_t         tofLowerQ;
  Float_t         tofQlay[4];

  int             nRichRing;
  float           richBeta;
  float           richBetaError;
  float           richProb;
  bool            richGoodBeta;
  bool            richIsGood;
  bool            richIsClean;
  float           richMass;
  bool            richIsNaF;
  float           richExpectedPhotoElectrons_true;
  float           richExpectedPhotoElectrons_false;
  float           richPhotoElectrons_true;
  float           richPhotoElectrons_false;
 
  Float_t         rtiLifetime;
  Float_t         rtiMaxCutoff[4][2];
  Int_t           passRTI;

  double BDT;
  double BDTG;
  
  //----------------------------------------------------

  // --- Register th trees
  TChain* TreeS   = new TChain("TreeS");  
  TChain* TreeST  = new TChain("TreeST");
  TChain* TreeB   = new TChain("TreeB");
  TChain* TreeBT  = new TChain("TreeBT");
  TChain* TreeSTD = new TChain("TreeSTD");
  TChain* TreeBTD = new TChain("TreeBTD");
  
  // Read training and test data
  printf("*** Reading input files from %s/\n", TrainingAndTestTreesPath.Data());
  for (int ii=bin_start; ii<=bin_stop; ii++) {
    TFile* input = new TFile(Form("%s/TrainingAndTestTrees_%02d_v%d.root", TrainingAndTestTreesPath.Data(), ii, SCversion));
    std::cout << "--- Plot       : Using input file: " << input->GetName() << std::endl;   
    if (!input) {
      cout << "ERROR: could not open data file" << endl;
    }
    else {
      TreeS->Add(input->GetName());
      TreeST->Add(input->GetName());
      TreeB->Add(input->GetName());
      TreeBT->Add(input->GetName());
      TreeSTD->Add(input->GetName());
      TreeBTD->Add(input->GetName());
      input->Close();
      delete input;
    }
  }
  
  //  TChain* ch[6] = {TreeS, TreeST, TreeB, TreeBT, TreeSTD, TreeBTD};
  TChain* ch[4] = {TreeS, TreeST, TreeB, TreeBT};
  Long64_t total_entries = 0; 
  //  Long64_t nentries[6] = {0, 0, 0, 0, 0, 0};
  Long64_t nentries[4] = {0, 0, 0, 0};
  
  for (int dd=0; dd<4; dd++) {
    printf("----------------------------------\n");
    printf("Chain %s\n", ch[dd]->GetName());
    printf("----------------------------------\n");
    
    ch[dd]->GetListOfFiles()->Print();
    
    nentries[dd] = ch[dd]->GetEntries();
    total_entries += nentries[dd];
    printf("TChain %s with %lld entries (%lld in total)\n", ch[dd]->GetName(), nentries[dd], total_entries);
    //  sleep(5);
  }

  unsigned int all_events = 0;
  unsigned int analyzed_events = 0;
  unsigned int stdcuts_events = 0;
  unsigned int addcut_events = 0;
  unsigned int signal_events = 0;

  for (int dd=0; dd<4; dd++) {
    
    ch[dd]->SetBranchStatus("*",0);  // disable all branches
    
    // Set branch addresses.
    SetBranchAddressAndActivate(ch[dd], "Run",           &Run);
    SetBranchAddressAndActivate(ch[dd], "Event",         &Event);
    SetBranchAddressAndActivate(ch[dd], "Time",         &Time);
    SetBranchAddressAndActivate(ch[dd], "UTime",        &UTime);
    SetBranchAddressAndActivate(ch[dd], "runAnalysisTag", &runAnalysisTag);
    
    SetBranchAddressAndActivate(ch[dd], "trdKlhrEP", &trdKlhrEP);
    SetBranchAddressAndActivate(ch[dd], "trdKlhrHP", &trdKlhrHP);
    SetBranchAddressAndActivate(ch[dd], "trdKlhrEH", &trdKlhrEH);
    // SetBranchAddressAndActivate(ch[dd], "trdKlhrEP_ene", &trdKlhrEP_ene);
    // SetBranchAddressAndActivate(ch[dd], "trdKlhrHP_ene", &trdKlhrHP_ene);
    // SetBranchAddressAndActivate(ch[dd], "trdKlhrEH_ene", &trdKlhrEH_ene);
    // SetBranchAddressAndActivate(ch[dd], "trdKlhE_ene", &trdKlhE_ene);
    // SetBranchAddressAndActivate(ch[dd], "trdKlhP_ene", &trdKlhP_ene);
    // SetBranchAddressAndActivate(ch[dd], "trdKlhH_ene", &trdKlhH_ene);
    // SetBranchAddressAndActivate(ch[dd], "trdKlhrEP_ene_trdstd", &trdKlhrEP_ene_trdstd);
    // SetBranchAddressAndActivate(ch[dd], "trdKlhrHP_ene_trdstd", &trdKlhrHP_ene_trdstd);
    // SetBranchAddressAndActivate(ch[dd], "trdKlhrEH_ene_trdstd", &trdKlhrEH_ene_trdstd);
    // SetBranchAddressAndActivate(ch[dd], "trdKlhE_ene_trdstd", &trdKlhE_ene_trdstd);
    // SetBranchAddressAndActivate(ch[dd], "trdKlhP_ene_trdstd", &trdKlhP_ene_trdstd);
    // SetBranchAddressAndActivate(ch[dd], "trdKlhH_ene_trdstd", &trdKlhH_ene_trdstd);
    SetBranchAddressAndActivate(ch[dd], "trdKHitsOnTrkTrack", &trdKHitsOnTrkTrack);
    SetBranchAddressAndActivate(ch[dd], "trdKHitsOnTrdTrack", &trdKHitsOnTrdTrack);
    
    SetBranchAddressAndActivate(ch[dd], "ecalBDTv5Es", &ecalBDTv5Es);
    SetBranchAddressAndActivate(ch[dd], "ecalBDTv6Es", &ecalBDTv6Es);
    SetBranchAddressAndActivate(ch[dd], "ecalBDTv7Es", &ecalBDTv7Es);
    SetBranchAddressAndActivate(ch[dd], "ecalCorrEneE_E",&ecalCorrEneE_E);
    SetBranchAddressAndActivate(ch[dd], "ecalCorrEneE_P", &ecalCorrEneE_P);
    SetBranchAddressAndActivate(ch[dd], "ecalCorrEneP_E", &ecalCorrEneP_E);
    SetBranchAddressAndActivate(ch[dd], "ecalCorrEneP_P", &ecalCorrEneP_P);
    
    SetBranchAddressAndActivate(ch[dd], "nTrTracks", &nTrTracks);
    SetBranchAddressAndActivate(ch[dd], "trkDefRig", &trkDefRig);
    SetBranchAddressAndActivate(ch[dd], "trkDefTheta", &trkDefTheta);
    SetBranchAddressAndActivate(ch[dd], "trkDefPhi", &trkDefPhi);
    SetBranchAddressAndActivate(ch[dd], "trkInnQ", &trkInnQ);
    
    SetBranchAddressAndActivate(ch[dd], "trigPhysBPatt", &trigPhysBPatt);

    SetBranchAddressAndActivate(ch[dd], "tofBeta", &tofBeta);
    SetBranchAddressAndActivate(ch[dd], "tofQ", &tofQ);
    SetBranchAddressAndActivate(ch[dd], "tofUpperQ", &tofUpperQ);
    SetBranchAddressAndActivate(ch[dd], "tofLowerQ", &tofLowerQ);
    SetBranchAddressAndActivate(ch[dd], "tofQlay", tofQlay);

    SetBranchAddressAndActivate(ch[dd], "nRichRing",       &nRichRing);
    SetBranchAddressAndActivate(ch[dd], "richBeta",        &richBeta);
    SetBranchAddressAndActivate(ch[dd], "richBetaError",   &richBetaError);
    SetBranchAddressAndActivate(ch[dd], "richProb",        &richProb);
    SetBranchAddressAndActivate(ch[dd], "richGoodBeta",    &richGoodBeta);
    SetBranchAddressAndActivate(ch[dd], "richIsGood",      &richIsGood);
    SetBranchAddressAndActivate(ch[dd], "richIsClean",     &richIsClean);
    SetBranchAddressAndActivate(ch[dd], "richMass",        &richMass);
    SetBranchAddressAndActivate(ch[dd], "richIsNaF",       &richIsNaF);
    SetBranchAddressAndActivate(ch[dd], "richExpectedPhotoElectrons_true",        &richExpectedPhotoElectrons_true);
    SetBranchAddressAndActivate(ch[dd], "richExpectedPhotoElectrons_false",       &richExpectedPhotoElectrons_false);
    SetBranchAddressAndActivate(ch[dd], "richPhotoElectrons_true",                &richPhotoElectrons_true);
    SetBranchAddressAndActivate(ch[dd], "richPhotoElectrons_false",               &richPhotoElectrons_false);
    
    //    if (dd==0) { //Data only
    SetBranchAddressAndActivate(ch[dd],  "rtiLifetime",&rtiLifetime);
    SetBranchAddressAndActivate(ch[dd],  "rtiMaxCutoff",rtiMaxCutoff);
    SetBranchAddressAndActivate(ch[dd],  "passRTI",&passRTI);
    //    }
    
    SetBranchAddressAndActivate(ch[dd],  "BDT", &BDT);
    SetBranchAddressAndActivate(ch[dd],  "BDTG", &BDTG);
    
    //----------------------------------------------------------------
    printf("All branches set\n");
    //----------------------------------------------------------------
      
    int perc=1;
  
    Long64_t nbytes = 0;
    
    bool kQUICK=false;

    double kPrescale = 1.0;
    double kAdditionalProtonPrescale = 1.0;
    
    Long64_t ee=0;
    for (ee=0; ee<nentries[dd]; ee++) {

      all_events++;
      
      if (gRandom->Uniform(0,1)>kPrescale) continue;
      
      nbytes += ch[dd]->GetEntry(ee);
      double pperc=(100.0*ee)/nentries[dd];
      if (pperc>perc) {
	printf("Processed %lld out of %lld: %d%%\n", ee+1, nentries[dd], (int)perc);
	perc++;
      }
      //      printf("Processing %lld...\n", ee);
      //      printf("Event=%u, Run=%u, Time=%u, UTime=%u\n", Event, Run, Time, UTime);
      //-------------------------------------------------------------------------------------------------------------------

      if (kQUICK && ee>100 ) break; 
    
      //-------------------------------------------------------------------------------------------------------------------
    
      bool cFisici = (trigPhysBPatt & 0x3E)!=0;
      bool cTrdHits = trdKHitsOnTrkTrack>8;
      bool cTrdDefined = trdKlhrEP>0;
      //      bool cTrdDefined = trdKlhE_ene>0;
      bool cCharge = trkInnQ<1.5;
      bool cTrdeHe = trdKlhrEH<0.8;
      //      bool cTrdeHe = trdKlhrEH_ene<0.8;
      bool cSingleTrack = nTrTracks<2;
      bool cGood = runAnalysisTag>=0 && passRTI>0 && rtiLifetime>0;
      //      bool cCutOff = ecalCorrEneE_E>1.25*TMath::Max(fabs(rtiMaxCutoff[0][0]), fabs(rtiMaxCutoff[0][1]));
      bool cCutOff = fabs(trkDefRig)>1.25*TMath::Max(fabs(rtiMaxCutoff[0][0]), fabs(rtiMaxCutoff[0][1]));
    
      bool cTrue = 1;
    
      bool stdcuts =
	cFisici &&
	cTrdHits &&
	cTrdDefined &&
	cCharge &&
	cTrdeHe &&
	cSingleTrack && 
	cGood &&
	cCutOff &&
	cTrue;
    
      //-------------------------------------------------------------------------------------------------------------------
    
      bool cNegative = trkDefRig<0.0;
      bool cPositive = trkDefRig>0.0;

      //-------------------------------------------------------------------------------------------------------------------

      bool cAddCut = fabs(trkDefRig)<25.0;

      bool cRICHCut = nRichRing>0 && richBeta>0;
      bool cECALCut = ecalCorrEneE_E>0.5;
    
      //      float trdclass = -(log10(trdKlhE_ene)+2);

      bool cElectronTRD = trdKlhrEP<0.5;
      bool cProtonTRD =  trdKlhrEP>0.8;
      
      bool cElectronECAL = cECALCut && (ecalBDTv5Es>0.2 && fabs(ecalCorrEneE_E/trkDefRig)>0.5);
      bool cElectronECALStrong = cECALCut && (ecalBDTv5Es>0.4 && fabs(ecalCorrEneE_E/trkDefRig)>0.75);
      bool cElectronRICH = cRICHCut && (BDT>0.4 && BDTG>0.3);
      bool cElectronRICHStrong = cRICHCut && (BDT>0.6 && BDTG>0.9);

      bool cProtonECAL = cECALCut && (ecalBDTv5Es<-0.4 && fabs(ecalCorrEneE_E/trkDefRig)<0.5);
      bool cProtonECALStrong = cECALCut && (ecalBDTv5Es<-0.8 && fabs(ecalCorrEneE_E/trkDefRig)<0.4);
      bool cProtonECALVeryStrong = cECALCut && (ecalBDTv5Es<-0.9 && fabs(ecalCorrEneE_E/trkDefRig)<0.3);
      bool cProtonECALHS = cECALCut && (ecalBDTv5Es<-0.5 && fabs(ecalCorrEneE_E/trkDefRig)<0.5);
      bool cProtonRICH = cRICHCut && (BDT<0.4 && BDTG<-0.2);
      bool cProtonRICHStrong = cRICHCut && (BDT<0.2 && BDTG<-0.8);
      bool cProtonRICHVeryStrong = cRICHCut && (BDT<0.0 && BDTG<-0.97);
      bool cProtonRICHHS = cRICHCut && (BDT<0.3 && BDTG<-0.5);

      //--------------------------------------------------------------------------------
      // bool cElectron = (trdKlhrEP_ene<0.5 && trkDefRig<0);
      // bool cProton =  (trdKlhrEP_ene>0.9 && trkDefRig>0); // ceci: trd killer cut
      //--------------------------------------------------------------------------------
      // bool cElectron = (trdKlhrEP_ene<0.6 && trkDefRig<0);
      // bool cProton =  (trdKlhrEP_ene>0.75 && trkDefRig>0);
      //--------------------------------------------------------------------------------
      // bool cElectron = (trdclass>1.0 && cNegative);
      // bool cProton = (trdclass<1.0 && cPositive);
      //--------------------------------------------------------------------------------
      //bool cElectron = (trdKlhrEP<0.5 && cNegative);
      //bool cProton =  (trdKlhrEP>0.8 && cPositive);
      //--------------------------------------------------------------------------------        
      bool cElectron = ((cElectronECAL || cElectronRICH) && cNegative);
      bool cElectronStrong = ((cElectronECALStrong || cElectronRICHStrong) && cNegative);
      bool cProton =  ((cProtonECAL || cProtonRICH) && cPositive);
      bool cProtonStrong =  ((cProtonECALStrong || cProtonRICHStrong) && cPositive);
      bool cProtonCC =  ((cProtonECAL || cProtonRICH) && cNegative);
      bool cProtonCCStrong =  ((cProtonECALStrong || cProtonRICHStrong) && cNegative);
      bool cProtonCCVeryStrong =  ((cProtonECALVeryStrong || cProtonRICHVeryStrong) && cNegative);
      bool cProtonCCHS =  ((cProtonECALHS || cProtonRICHHS) && cNegative);
      //--------------------------------------------------------------------------------        

      //      printf("%d %d %d\n", cElectron, cProton, cProtonCC);
      
      // printf("Energy = %f\n", ecalCorrEneE_E);

      analyzed_events++;
      
      if (stdcuts) {

	stdcuts_events++;
	
	if (cAddCut) { 
	  //	  printf("Event=%u, Run=%u: passed std-cuts\n", Event, Run);
	  addcut_events++;
	  
	  //---------------------------------------------------------------
		
	  hEventsAll->Fill(Time);

	  if (cRICHCut) {
	    if (cNegative) {
	      hBDTvsRNegative->Fill(log10(fabs(trkDefRig)), BDT);
	      hBDTGvsRNegative->Fill(log10(fabs(trkDefRig)), BDTG);
	      if (cProtonTRD) {
		hBDTvsRProtonsCC->Fill(log10(fabs(trkDefRig)), BDT);
		hBDTGvsRProtonsCC->Fill(log10(fabs(trkDefRig)), BDTG);
	      }
	      if (cElectronTRD) {
		hBDTvsRElectrons->Fill(log10(fabs(trkDefRig)), BDT);
		hBDTGvsRElectrons->Fill(log10(fabs(trkDefRig)), BDTG);
	      }
	    }
	    if (cPositive) {
	      hBDTvsRPositive->Fill(log10(fabs(trkDefRig)), BDT);
	      hBDTGvsRPositive->Fill(log10(fabs(trkDefRig)), BDTG);
	      if (cProtonTRD) {
		hBDTvsRProtons->Fill(log10(fabs(trkDefRig)), BDT);
		hBDTGvsRProtons->Fill(log10(fabs(trkDefRig)), BDTG);
	      }
	    }
	  }
	    
	  if (!(cElectron || cProton || cProtonCC)) continue;
	  if (cProton && gRandom->Uniform(0,1)>kAdditionalProtonPrescale) continue;

	  signal_events++;
	  
	  //---------------------------------------------------------------

	  double mu=0;
	  double sigma=1;
	  double targetmu = 1.05;
	  double targetsigma = 0.175;

	  if (_trdmap) {
	    std::pair<double, double> meansigma = _trdmap->getelementbytime(Time);
	    //		printf("Time = %u, Mean=%f\n", Time, meansigma.first);
	    mu = meansigma.first;
	    sigma = meansigma.second;
	  }
	  
	  hEventsRICHorECAL->Fill(Time);
	  hTRDAll->Fill(trdKlhrEP);
	  if (cECALCut) {
	    hEventsECAL->Fill(Time);
	    hTRDECAL->Fill(trdKlhrEP);
	  }
	  if (cRICHCut) {
	    hEventsRICH->Fill(Time);
	    hTRDRICH->Fill(trdKlhrEP);
	  }

	  if (cElectron) {
	    hEventsElectrons->Fill(Time);
	    if (cElectronStrong) {
	      hTRDElectrons->Fill(trdKlhrEP);
	      hTRDElectronsVsTime->Fill(Time, trdKlhrEP);
	      hTRDNormalizedElectronsVsTime->Fill(Time, targetmu+targetsigma*(trdKlhrEP-mu)/sigma);
	      hTRDvsRElectrons->Fill(log10(fabs(trkDefRig)), trdKlhrEP);
	    }
	    if (cECALCut) {
	      hEventsElectronsECAL->Fill(Time);
	      if (cElectronECALStrong && cNegative) {
		hTRDElectronsECAL->Fill(trdKlhrEP);
		hTRDvsRElectronsECAL->Fill(log10(fabs(trkDefRig)), trdKlhrEP);
	      }
	    }
	    if (cRICHCut) {
	      hEventsElectronsRICH->Fill(Time);
	      if (cElectronRICHStrong && cNegative) {
		hTRDElectronsRICH->Fill(trdKlhrEP);
		hTRDvsRElectronsRICH->Fill(log10(fabs(trkDefRig)), trdKlhrEP);
	      }
	    }
	  }
	  else if (cProton) {
	    hEventsProtons->Fill(Time);
	    if (cProtonStrong) {
	      hTRDProtons->Fill(trdKlhrEP);
	      hTRDProtonsVsTime->Fill(Time, trdKlhrEP);
	      hTRDNormalizedProtonsVsTime->Fill(Time, targetmu+targetsigma*(trdKlhrEP-mu)/sigma);
	      hTRDvsRProtons->Fill(log10(fabs(trkDefRig)), trdKlhrEP);
	    }
	    if (cECALCut) {
	      hEventsProtonsECAL->Fill(Time);
	      if (cProtonECALStrong && cPositive) {
		hTRDProtonsECAL->Fill(trdKlhrEP);
		hTRDvsRProtonsECAL->Fill(log10(fabs(trkDefRig)), trdKlhrEP);
	      }
	    }
	    if (cRICHCut) {
	      hEventsProtonsRICH->Fill(Time);
	      if (cProtonRICHStrong && cPositive) {
		hTRDProtonsRICH->Fill(trdKlhrEP);
		hTRDvsRProtonsRICH->Fill(log10(fabs(trkDefRig)), trdKlhrEP);
	      }
	    }
	  }
	  else if (cProtonCC) {
	    if (cProtonCCHS) {
	      hTRDProtonsCCHS->Fill(trdKlhrEP);
	      hTRDProtonsCCHSVsTime->Fill(Time, trdKlhrEP);
	      hTRDNormalizedProtonsCCHSVsTime->Fill(Time, targetmu+targetsigma*(trdKlhrEP-mu)/sigma);
	      hTRDvsRProtonsCCHS->Fill(log10(fabs(trkDefRig)), trdKlhrEP);
	    }
	    if (cProtonCCVeryStrong) {
	      hTRDProtonsCC->Fill(trdKlhrEP);
	      hTRDProtonsCCVsTime->Fill(Time, trdKlhrEP);
	      hTRDNormalizedProtonsCCVsTime->Fill(Time, targetmu+targetsigma*(trdKlhrEP-mu)/sigma);
	      hTRDvsRProtonsCC->Fill(log10(fabs(trkDefRig)), trdKlhrEP);
	    }
	    if (cECALCut) {
	      if (cProtonECALHS && cNegative) {
		hTRDProtonsCCHSECAL->Fill(trdKlhrEP);
		hTRDvsRProtonsCCHSECAL->Fill(log10(fabs(trkDefRig)), trdKlhrEP);
	      }
	      if (cProtonECALVeryStrong && cNegative) {
		hTRDProtonsCCECAL->Fill(trdKlhrEP);
		hTRDvsRProtonsCCECAL->Fill(log10(fabs(trkDefRig)), trdKlhrEP);
	      }
	    }
	    if (cRICHCut) {
	      if (cProtonRICHHS && cNegative) {
		hTRDProtonsCCHSRICH->Fill(trdKlhrEP);
		hTRDvsRProtonsCCHSRICH->Fill(log10(fabs(trkDefRig)), trdKlhrEP);
	      }
	      if (cProtonRICHVeryStrong && cNegative) {
		hTRDProtonsCCRICH->Fill(trdKlhrEP);
		hTRDvsRProtonsCCRICH->Fill(log10(fabs(trkDefRig)), trdKlhrEP);
	      }
	    }
	  }
	  
	  if (cElectron) {
	    if (fabs(180.0*trkDefTheta/TMath::Pi())<5.0)       hTRDElectrons0005->Fill(trdKlhrEP);
	    else if (fabs(180.0*trkDefTheta/TMath::Pi())<10.0) hTRDElectrons0510->Fill(trdKlhrEP);
	    else if (fabs(180.0*trkDefTheta/TMath::Pi())<15.0) hTRDElectrons1015->Fill(trdKlhrEP);
	    else if (fabs(180.0*trkDefTheta/TMath::Pi())<20.0) hTRDElectrons1520->Fill(trdKlhrEP);
	    else if (fabs(180.0*trkDefTheta/TMath::Pi())<25.0) hTRDElectrons2025->Fill(trdKlhrEP);
	    else if (fabs(180.0*trkDefTheta/TMath::Pi())<30.0) hTRDElectrons2530->Fill(trdKlhrEP);
	    else if (fabs(180.0*trkDefTheta/TMath::Pi())<35.0) hTRDElectrons3035->Fill(trdKlhrEP);
	  }

	  if (cProton) {
	    if (fabs(180.0*trkDefTheta/TMath::Pi())<5.0)       hTRDProtons0005->Fill(trdKlhrEP);
	    else if (fabs(180.0*trkDefTheta/TMath::Pi())<10.0) hTRDProtons0510->Fill(trdKlhrEP);
	    else if (fabs(180.0*trkDefTheta/TMath::Pi())<15.0) hTRDProtons1015->Fill(trdKlhrEP);
	    else if (fabs(180.0*trkDefTheta/TMath::Pi())<20.0) hTRDProtons1520->Fill(trdKlhrEP);
	    else if (fabs(180.0*trkDefTheta/TMath::Pi())<25.0) hTRDProtons2025->Fill(trdKlhrEP);
	    else if (fabs(180.0*trkDefTheta/TMath::Pi())<30.0) hTRDProtons2530->Fill(trdKlhrEP);
	    else if (fabs(180.0*trkDefTheta/TMath::Pi())<35.0) hTRDProtons3035->Fill(trdKlhrEP);
	  }

	}
      }
    
    }//loop on entries
    printf("All entries looped...\n");
  
  }//loop on chains
  printf("All input chains processed...\n");

  //---------------------------------------------------
  
  printf("all events: %u\n", all_events);
  printf("analyzed events: %u\n", analyzed_events);
  printf("stdcuts events: %u\n", stdcuts_events);
  printf("addcut events: %u\n", addcut_events);
  printf("signal events: %u\n", signal_events);
  
  //---------------------------------------------------
  
  // Save the output
  outputFile->Write(outputFile->GetName(), TObject::kOverwrite);  
  std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
  outputFile->Close();

  return;
}

// -------------------------------------------------------------------------------------

TH2* SliceNormalizeX(TH2* h){
  
  if (!h) return 0;
  
  TH2* hclone = (TH2*)(h->Clone(Form("%s_SliceNormalizedX", h->GetName())));

  for (int ii=0; ii<=h->GetNbinsX()+1; ii++) {
    //    printf("Binx %d\n", ii);
    double integralslice=0;
    integralslice=h->Integral(ii, ii, 0, h->GetNbinsY()+1);
    //    printf("Integral = %f\n", integralslice);
    if (integralslice>0) {
      for (int jj=0; jj<=h->GetNbinsY()+1; jj++) {
        if (integralslice) {
          hclone->SetBinContent(ii, jj, h->GetBinContent(ii, jj)/integralslice);
          hclone->SetBinError(ii, jj, h->GetBinError(ii, jj)/integralslice);
        }
      }
    }
  }

  return hclone;
}

TH2* SliceNormalizeY(TH2* h){
  
  if (!h) return 0;

  TH2* hclone = (TH2*)(h->Clone(Form("%s_SliceNormalizedY", h->GetName())));  

  for (int ii=0; ii<=h->GetNbinsY()+1; ii++) {
    //    printf("Biny %d\n", ii);
    double integralslice=0;
    integralslice=h->Integral(0, h->GetNbinsX()+1, ii, ii);
    //    printf("Integral = %f\n", integralslice);
    if (integralslice>0) {
      for (int jj=0; jj<=h->GetNbinsX()+1; jj++) {
        if (integralslice) {
          hclone->SetBinContent(jj, ii, h->GetBinContent(jj, ii)/integralslice);
          hclone->SetBinError(jj, ii, h->GetBinError(jj, ii)/integralslice);
        }
      }
    }
  }

  return hclone;
}

TH1* CreateHistoFromGraph(TGraph* gr){

  if (!gr) return 0;

  TH1F* hh = new TH1F(Form("%s_histo", gr->GetName()), gr->GetTitle(), nbozzo-1, bozzo_bins);
  
  Double_t X,Y;
  
  for(int i=0; i<gr->GetN(); i++) {
    gr->GetPoint(i,X,Y);
    hh->SetBinContent(hh->GetXaxis()->FindBin(X), Y);
  }

  return hh;
}
