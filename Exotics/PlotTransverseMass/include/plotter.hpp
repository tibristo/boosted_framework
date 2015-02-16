#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TGraph.h"
#include "Rtypes.h"
#include <string>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <map>
#include <cmath>
#include "TH1F.h"
#include "TROOT.h"
#include "TInterpreter.h"
#include "TPostScript.h"
#include "TCanvas.h"
#include "TTreeIndex.h"
#include "TMath.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TLorentzVector.h"


class plotter//;
{
public:
  plotter(std::string & algorithm, TTree * bkg_in, TTree * sig_in, std::string & alg_params, std::string & dir, TH1F * sig_hist, TH1F * bkg_hist, bool applyPt, bool applyEta);
  void returnROC();
  void createHistograms(std::map<std::string, std::vector<std::string> > & config);
  void getVariables();
  void fillHistograms(TTree * tree, bool signal);
  void plotHistograms();
  void readAlgorithmConfig(std::string &algorithm);
  TGraph produceROC(int type, TH1F *&S, TH1F *&B, TString name, bool draw);
  void initVectors();
  void setBranches();
  void runPlotter();
  void clearHists();
  void initHists();
  void produceAllROCs(std::string alg);
  void getWeights(std::string & algorithm, std::string & params);
  void compareHistograms(TH1F * hist1, TH1F * hist2, std::string file_name);

  float getPtWeight(double & pt);

  std::string algorithm;
  std::string algorithm_params;
  std::string working_dir;

  enum jet_type {TRUTH, TOPO, GROOMED};

  TTree * sig;
  TTree * bkg;

  bool applyPtCut;
  bool applyEtaCut;

  TGraph tg;

  Int_t signal_truth_index;
  Int_t signal_topo_index;
  Int_t signal_groomed_index;

  Int_t bkg_truth_index;
  Int_t bkg_topo_index;
  Int_t bkg_groomed_index;

  Float_t sig_normalisation = 1.0;
  Float_t bkg_normalisation = 1.0;
  Float_t sig_ptreweight = 1.0;
  Float_t bkg_ptreweight = 1.0;
  Float_t signal_mc_event_number = 0;
  Float_t bkg_mc_event_number = 0;
  Float_t signal_mc_event_weight = 1.0;
  Float_t bkg_mc_event_weight = 1.0;
  //Float_t signal_NEvents_weighted = 1.0;
  //Float_t bkg_NEvents_weighted = 1.0;
  Int_t signal_NEvents = 1.0;
  Int_t bkg_NEvents = 1.0;
  UInt_t signal_RunNumber = 0;
  UInt_t bkg_RunNumber = 0;
  std::map<int,float> events;


  TH1F * Wp_PtReweight;
  TH1F * qcd_PtReweight;
  TH1F * pt_reweight_hist;

  std::map<std::string, std::vector<std::string> > config;
  std::map<Int_t, std::vector<std::string> > weighting;
  std::map<std::string, std::string> vars;

  std::vector<TH1F *> signal_E;
  std::vector<TH1F *> signal_pt;
  std::vector<TH1F *> signal_m;
  std::vector<TH1F *> signal_eta;
  std::vector<TH1F *> signal_phi;
  std::vector<TH1F *> signal_emfrac;
  std::vector<TH1F *> signal_Tau1;
  std::vector<TH1F *> signal_Tau2;
  std::vector<TH1F *> signal_Tau3;
  std::vector<TH1F *> signal_WIDTH;
  std::vector<TH1F *> signal_SPLIT12;
  std::vector<TH1F *> signal_SPLIT23;
  std::vector<TH1F *> signal_SPLIT34;
  std::vector<TH1F *> signal_Dip12;
  std::vector<TH1F *> signal_Dip13;
  std::vector<TH1F *> signal_Dip23;
  std::vector<TH1F *> signal_DipExcl12;
  std::vector<TH1F *> signal_PlanarFlow;
  std::vector<TH1F *> signal_Angularity;
  std::vector<TH1F *> signal_QW;
  std::vector<TH1F *> signal_PullMag;
  std::vector<TH1F *> signal_PullPhi;
  std::vector<TH1F *> signal_Pull_C00;
  std::vector<TH1F *> signal_Pull_C01;
  std::vector<TH1F *> signal_Pull_C10;
  std::vector<TH1F *> signal_Pull_C11;


  std::vector<TH1F *> bkg_E;
  std::vector<TH1F *> bkg_pt;
  std::vector<TH1F *> bkg_m;
  std::vector<TH1F *> bkg_eta;
  std::vector<TH1F *> bkg_phi;
  std::vector<TH1F *> bkg_emfrac;
  std::vector<TH1F *> bkg_Tau1;
  std::vector<TH1F *> bkg_Tau2;
  std::vector<TH1F *> bkg_Tau3;
  std::vector<TH1F *> bkg_WIDTH;
  std::vector<TH1F *> bkg_SPLIT12;
  std::vector<TH1F *> bkg_SPLIT23;
  std::vector<TH1F *> bkg_SPLIT34;
  std::vector<TH1F *> bkg_Dip12;
  std::vector<TH1F *> bkg_Dip13;
  std::vector<TH1F *> bkg_Dip23;
  std::vector<TH1F *> bkg_DipExcl12;
  std::vector<TH1F *> bkg_PlanarFlow;
  std::vector<TH1F *> bkg_Angularity;
  std::vector<TH1F *> bkg_QW;
  std::vector<TH1F *> bkg_PullMag;
  std::vector<TH1F *> bkg_PullPhi;
  std::vector<TH1F *> bkg_Pull_C00;
  std::vector<TH1F *> bkg_Pull_C01;
  std::vector<TH1F *> bkg_Pull_C10;
  std::vector<TH1F *> bkg_Pull_C11;


  // have the vectors for the above histograms so we can do the reading in stuff from the TTree
  std::map<int, std::vector<Float_t> *> signal_E_vec;
  std::map<int, std::vector<Float_t> *> signal_pt_vec;
  std::map<int, std::vector<Float_t> *> signal_m_vec;
  std::map<int, std::vector<Float_t> *> signal_eta_vec;
  std::map<int, std::vector<Float_t> *> signal_phi_vec;
  std::map<int, std::vector<Float_t> *> signal_emfrac_vec;
  std::map<int, std::vector<Float_t> *> signal_Tau1_vec;
  std::map<int, std::vector<Float_t> *> signal_Tau2_vec;
  std::map<int, std::vector<Float_t> *> signal_Tau3_vec;
  std::map<int, std::vector<Float_t> *> signal_WIDTH_vec;
  std::map<int, std::vector<Float_t> *> signal_SPLIT12_vec;
  std::map<int, std::vector<Float_t> *> signal_SPLIT23_vec;
  std::map<int, std::vector<Float_t> *> signal_SPLIT34_vec;
  std::map<int, std::vector<Float_t> *> signal_Dip12_vec;
  std::map<int, std::vector<Float_t> *> signal_Dip13_vec;
  std::map<int, std::vector<Float_t> *> signal_Dip23_vec;
  std::map<int, std::vector<Float_t> *> signal_DipExcl12_vec;
  std::map<int, std::vector<Float_t> *> signal_PlanarFlow_vec;
  std::map<int, std::vector<Float_t> *> signal_Angularity_vec;
  std::map<int, std::vector<Float_t> *> signal_QW_vec;
  std::map<int, std::vector<Float_t> *> signal_PullMag_vec;
  std::map<int, std::vector<Float_t> *> signal_PullPhi_vec;
  std::map<int, std::vector<Float_t> *> signal_Pull_C00_vec;
  std::map<int, std::vector<Float_t> *> signal_Pull_C01_vec;
  std::map<int, std::vector<Float_t> *> signal_Pull_C10_vec;
  std::map<int, std::vector<Float_t> *> signal_Pull_C11_vec;



  std::map<int, std::vector<Float_t> *> bkg_E_vec;
  std::map<int, std::vector<Float_t> *> bkg_pt_vec;
  std::map<int, std::vector<Float_t> *> bkg_m_vec;
  std::map<int, std::vector<Float_t> *> bkg_eta_vec;
  std::map<int, std::vector<Float_t> *> bkg_phi_vec;
  std::map<int, std::vector<Float_t> *> bkg_emfrac_vec;
  std::map<int, std::vector<Float_t> *> bkg_Tau1_vec;
  std::map<int, std::vector<Float_t> *> bkg_Tau2_vec;
  std::map<int, std::vector<Float_t> *> bkg_Tau3_vec;
  std::map<int, std::vector<Float_t> *> bkg_WIDTH_vec;
  std::map<int, std::vector<Float_t> *> bkg_SPLIT12_vec;
  std::map<int, std::vector<Float_t> *> bkg_SPLIT23_vec;
  std::map<int, std::vector<Float_t> *> bkg_SPLIT34_vec;
  std::map<int, std::vector<Float_t> *> bkg_Dip12_vec;
  std::map<int, std::vector<Float_t> *> bkg_Dip13_vec;
  std::map<int, std::vector<Float_t> *> bkg_Dip23_vec;
  std::map<int, std::vector<Float_t> *> bkg_DipExcl12_vec;
  std::map<int, std::vector<Float_t> *> bkg_PlanarFlow_vec;
  std::map<int, std::vector<Float_t> *> bkg_Angularity_vec;
  std::map<int, std::vector<Float_t> *> bkg_QW_vec;
  std::map<int, std::vector<Float_t> *> bkg_PullMag_vec;
  std::map<int, std::vector<Float_t> *> bkg_PullPhi_vec;
  std::map<int, std::vector<Float_t> *> bkg_Pull_C00_vec;
  std::map<int, std::vector<Float_t> *> bkg_Pull_C01_vec;
  std::map<int, std::vector<Float_t> *> bkg_Pull_C10_vec;
  std::map<int, std::vector<Float_t> *> bkg_Pull_C11_vec;




  
};
static inline std::string &rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
  }
