
#include <TH1F.h>
#include "TChain.h"
#include <TH1I.h>
#include "TH1.h"
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

#include "TROOT.h"
#include "TInterpreter.h"
#include "TPostScript.h"
#include "TCanvas.h"
#include "TTreeIndex.h"
#include "TMath.h"
#include "AtlasStyle.C"
#include "TLegend.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TLorentzVector.h"
//FASTJET INCLUDES

#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequence1GhostPassiveArea.hh"
#include "fastjet/ClusterSequenceActiveArea.hh"
#include "fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/ClusterSequenceAreaBase.hh"
#include "fastjet/ClusterSequencePassiveArea.hh"
#include "fastjet/ClusterSequenceStructure.hh"
#include "fastjet/ClusterSequenceVoronoiArea.hh"
#include "fastjet/CompositeJetStructure.hh"
#include "fastjet/EECambridgePlugin.hh"
#include "fastjet/Error.hh"
#include "fastjet/FunctionOfPseudoJet.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/GridJetPlugin.hh"
#include "fastjet/JadePlugin.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/LimitedWarning.hh"
#include "fastjet/NNH.hh"
#include "fastjet/NestedDefsPlugin.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/PseudoJetStructureBase.hh"
#include "fastjet/RangeDefinition.hh"
#include "fastjet/SISConeBasePlugin.hh"
#include "fastjet/SISConePlugin.hh"
#include "fastjet/SISConeSphericalPlugin.hh"
#include "fastjet/Selector.hh"
#include "fastjet/SharedPtr.hh"
#include "fastjet/WrappedStructure.hh"
#include "fastjet/version.hh"



//Declarations
//bool CheckBranchExists( string BranchName );
float DeltaR(float eta1,float phi1,float eta2,float phi2);
double mpv(TH1F* histo);
void Qw(double &minWidth, double &topEdge, TH1F* histo, double frac);

//void cleanBranches();
//void defineStrings(TString *AlgoList, TString *binLabel, TString *pTbins, TString *finePtBins);
void createHistos(std::string & algo);
void getMassHistograms(TTree *inputTree, TTree *inputTree1, TString groomAlgo, std::string & groomAlgoIndex);
void initializeVariables();
void getBranches(TTree *inputTree, TTree *inputTree1, TString groomAlgo, std::string & groomAlgoIndex);
void deleteVectors();
void getNormSherpaW(TString inputTree, unsigned long & evnum, double & weight);
void makeROC(int type, TH1F *&S,TH1F *&B,TGraph &curve, TString name="", bool draw=false);
void makePlots(std::string & algorithm);
void makePtPlots(std::string & algorithm);
void setMassBranch(TTree * tree, std::string &algorithm, std::string & groomAlgoIndex);
void plotVariables(TTree * tree, vector<std::string> & branches);
std::vector<std::string> getListOfBranches(std::string &algorithm);
//void make68Plots(int algidx, TChain * bkg, TChain * sig);
void makeMassWindowFile(bool applyMassWindow, std::string & algorithm);
void setJetsBranches(TChain * tree, std::string & algorithm, std::string & groomIdx);
void addInfoBranches(TTree * tree);
void addSubJets(TTree * tree, std::string & algorithm, std::string & groomIdx);
//void getBranchesSelection(TTree * tree, std::string & algorithm);
void setSelectionVectors();
void runAlgorithm(TChain *inputTree, TChain *inputTree1, TString groomAlgo, std::string & groomAlgoIndex, bool massHistos);
vector<std::pair<std::string,int> > getListOfJetBranches(std::string & algorithm, TObjArray * brancharray);
void initVectors();
std::pair<int,int> getTwoLeadingSubjets(std::vector<int> & indices, std::vector<float> * subjets);
std::string returnJetType(std::string & samplePrefix, std::string & groomalgo, bool addLC, int idx);
std::string returnSubJetType(std::string & samplePrefix, std::string & groomalgo, bool addLC);
void clearVectors();
void clearOutputVariables();
void setOutputVariables(int idx1, int idx2, int idx3, int subidx, std::string & groomalgo, std::string & groomIdx);
void setOutputBranches(TTree* tree, std::string & algorithm, std::string & groomIdx);
void resetOutputVariables();
void getMPV();
void scaleHists();
void setAddress(TChain * tree, std::string  name, std::vector<Float_t> * var_vec);
void readWeights();
void createPtReweightFile(TH1F * bkg, TH1F * sig, std::string & fname);
void overlapRemoval();
int eventSelection();
bool leptonSelection(int lepType);
void setLeptons(TChain * tree, TObjArray * list);
void addLeptonBranches(std::string & jetString, TChain * tree);
void setLeptonVectors();
std::vector<float> dummyCharge(int size);
void setVector(TChain *& tree, TObjArray *& list, vector<TLorentzVector> *& vec, std::string branch);//, const char * branch);//std::string & branch);

void setVector(TChain *& tree, TObjArray *& list, vector<Float_t> *& vec, std::string branch);//, const char * branch);//std::string & branch);
bool useBranch(std::string branch, bool partialmatch = false);
void setLLJMass(int jetidx);
double calculateFoxWolfram20(vector<TLorentzVector>& clusters);
int calculateSoftDropTag(vector<TLorentzVector> & cluster);
vector<double> calculateQJetsVol(vector<TLorentzVector> & EventParticles, float radius);
std::vector<TLorentzVector> createClusters(int jettype, int jetidx);
double calculateEEC(int jettype, float beta=0.3, float exp=2);
void setRadius(std::string & prefix);


enum class groomAlgoEnum{groomZero, TopoSplitFilteredMu67SmallR0YCut9, TopoSplitFilteredMu100SmallR30YCut4, TopoTrimmedPtFrac5SmallR30, TopoTrimmedPtFrac5SmallR20, TopoPrunedCaRcutFactor50Zcut10, TopoPrunedCaRcutFactor50Zcut20, AntiKt2LCTopo, AntiKt3LCTopo, AntiKt4LCTopo};
enum sampleType{BACKGROUND,SIGNAL};
enum jetType{TRUTH,TOPO,GROOMED,MAX};
enum histType{TRUTHJET,GROOMEDJET,LEADTRUTHJET};
enum leptonType{FAIL,ELECTRON,MUON};

vector<int> algoMap; // stores the algorithm used for inputTree[x] in the case that we are not using all the input types - split/filter, trim, prune and recluster
std::map<int, int> fileMap; // stores the index of the input file for each grooming algorithm

//DECLARATIONS FOR RECLUSTERING FUNCTIONS

vector<fastjet::PseudoJet> ObjsToPJ(vector<TLorentzVector> jets);
vector<TLorentzVector> Recluster(vector<TLorentzVector> small_jets, double PTcut=15, double fcut=0.05, double jetRad=1.0);

struct Algorithms
{
  std::map<std::string, std::string> AlgoNames; // the jet groomieng alg
  std::map<std::string, std::string> AlgoPrefix; // the jet reco algorithm
  std::map<std::string, std::string> AlgoType; // if split/filter, trimmed, etc
  std::map<std::string, std::string> AlgoList; // the abbreviated name
  std::map<std::string, std::string> AlgoListN; // the plotting labels
  std::map<std::string, std::string> subjetMap; // subjet grooming alg for jet grooming alg 
  std::map<std::string, std::string> subjetIndex;
  std::map<std::string, std::string> binLabel;
  void load(const std::string & filename);
};

struct Algorithms algorithms;
/////

std::string fileid_global;
std::string treeName;
std::string branchesFile;
std::map<string, int> branchmap;
float GEV = 1000.;
float ELMASS = 0.511;
float MUMASS = 105.7;

bool calcQJets = false;
bool calcFoxWolfram20 = false;
bool calcSoftDrop = false;
bool calcEEC = false;

float radius = 1.0;

TFile *inputFile[2];
TTree *inputTree[2];
TChain *inputTChain[2];

//get the branches we want to use
Int_t leadGroomedIndex = 0;
Int_t leadTruthIndex = 0;
Int_t leadTopoIndex = 0;
TH1F * pt_reweight = 0;
TH1F * pt_reweight_arr[2];
Float_t normalisation = 1.0;
Int_t NEvents = 0;
std::map<int, float> NEvents_weighted;
Float_t mc_event_weight = 1.0;
Float_t mc_event_weight_out = 1.0;
UInt_t mc_channel_number  = 0;
UInt_t mc_channel_number_out  = 0;
UInt_t runNumberOut = 0;
UInt_t runNumberIn = 0;
Float_t avgIntpXingOut = 0;
Float_t avgIntpXingIn = 0;
// xaod expects UInt_t, but D3PDs Int_t. This is annoying and will cause problems
UInt_t nvtxIn = 0;
UInt_t nvtxOut = 0;
bool xAODJets = false;
bool xAODemfrac = false;
bool hvtllqq = false;
//QCD split filtering with Y cut 9
vector<float> * qcd_CA12_truth_pt = 0;
vector<float> * qcd_CA12_truth_eta = 0;
vector<float> * qcd_CA12_truth_phi = 0;
vector<float> * qcd_CA12_truth_mass = 0;
vector<float> * qcd_CA12_truth_E = 0;

vector<float> * qcd_CA12_topo_pt = 0;
vector<float> * qcd_CA12_topo_eta = 0;
vector<float> * qcd_CA12_topo_phi = 0;
vector<float> * qcd_CA12_topo_mass = 0;
vector<float> * qcd_CA12_topo_E = 0;
vector<float> * qcd_CA12_topo_emfrac = 0;

vector<float> * qcd_CA12_groomed_pt = 0;
vector<float> * qcd_CA12_groomed_eta = 0;
vector<float> * qcd_CA12_groomed_phi = 0;
vector<float> * qcd_CA12_groomed_mass = 0;
vector<float> * qcd_CA12_groomed_E = 0;
vector<float> * qcd_CA12_groomed_emfrac = 0;


vector<float> * jet_eta_truth = 0;
vector<float> * jet_phi_truth = 0;
//vector<float> * jet_emfrac_truth = 0;
vector<float> * jet_pt_truth = 0;
vector<float> * jet_m_truth = 0;

vector<float> * jet_eta_topo = 0;
vector<float> * jet_phi_topo = 0;
vector<float> * jet_emfrac_topo = 0;
vector<float> * jet_pt_topo = 0;

vector<float> * jet_eta_groomed = 0;
vector<float> * jet_phi_groomed = 0;
vector<float> * jet_emfrac_groomed = 0;
vector<float> * jet_pt_groomed = 0;


UInt_t qcd_mc_channel_number = 0; 
Float_t qcd_mc_event_weight = 0; 



//Wprime split filtering both cuts

vector<float> * Wp_CA12_truth_pt = 0;
vector<float> * Wp_CA12_truth_eta = 0;
vector<float> * Wp_CA12_truth_phi = 0;
vector<float> * Wp_CA12_truth_mass = 0;
vector<float> * Wp_CA12_truth_E = 0;
UInt_t Wp_mc_channel_number = 0; 
Float_t Wp_mc_event_weight = 0; 

vector<float> * Wp_CA12_topo_pt = 0;
vector<float> * Wp_CA12_topo_eta = 0;
vector<float> * Wp_CA12_topo_phi = 0;
vector<float> * Wp_CA12_topo_mass = 0;
vector<float> * Wp_CA12_topo_E = 0;
vector<float> * Wp_CA12_topo_emfrac = 0;

vector<float> * Wp_CA12_groomed_pt = 0;
vector<float> * Wp_CA12_groomed_eta = 0;
vector<float> * Wp_CA12_groomed_phi = 0;
vector<float> * Wp_CA12_groomed_mass = 0;
vector<float> * Wp_CA12_groomed_E = 0;
vector<float> * Wp_CA12_groomed_emfrac = 0;

vector<float> * currTree_mass = 0;

// create a hashtable/ map for all the variables we want to plot.... Index on variable name - we can read this in from xml config

//counters for efficiency issues
int nEvt_0=0; //total events
int nEvt_1=0; // events with leading proper topo jet
int nEvt_2=0; // events with matched truth to proper topo jet
int nEvt_3=0; // events with split cut 9 matched to truth jet matches to proper topo jet
int nEvt_4=0;

//counters for efficiency issues
int nEvt1_0=0; //total events
int nEvt1_1=0; // events with leading proper topo jet
int nEvt1_2=0; // events with matched truth to proper topo jet
int nEvt1_3=0; // events with split cut 9 matched to truth jet matches to proper topo jet
int nEvt1_4=0;


const int nAlgosMax=12;
int nAlgos=0;
//TString AlgoList[nAlgosMax];
//std::string AlgoNames[nAlgosMax];
//std::map<std::string, std::string> subjetMap;
//std::string AlgoPrefix[nAlgosMax];
//std::map<std::string, std::string> subjetIndex;
//TString binLabel[nAlgosMax-2];
const int nPtBins=6;
TString pTbins[nPtBins];
const int nFineBins=12;
TString finePtBins[nFineBins];
std::vector<pair<float,float> > ptrange;

TH1F *qcd_finePtBin_mass[3][nFineBins];
TH1F *Wprime_finePtBin_mass[3][nFineBins];

TH1F *qcd_Lead_CA12_mass[3][nPtBins]; 
TH1F *Wprime_Lead_CA12_mass[3][nPtBins];
TH1F *qcd_Lead_CA12_pt[3];
TH1F *Wprime_Lead_CA12_pt[3];
TH1F *Wprime_Lead_CA12_scaled_pt;
TH1F *pTweights;
TH1F *qcd_PtReweight;
TH1F *Wp_PtReweight;
TString AlgoListN;
TString pTbinsN[nPtBins];
TString pTFinebinsN[nFineBins];

double myMPV[nPtBins];
double WidthMassWindow[nPtBins];
double TopEdgeMassWindow[nPtBins];
double BottomEdgeMassWindow[nPtBins];
double QCDfrac[nPtBins];

double myMPV_finePt[nFineBins];
double WidthMassWindow_finePt[nFineBins];
double TopEdgeMassWindow_finePt[nFineBins];
double BottomEdgeMassWindow_finePt[nFineBins];
double QCDfrac_finePt[nFineBins];


TH1F * hMassLow[nPtBins];
TH1F * hMassHigh[nPtBins];
TH1F * hMPV[nPtBins];
TH1F * hWmass[nPtBins];
TH1F * hQCDeff[nPtBins];

TH1F * hQCDeff_finePt;


TH1F * windowsVsPt;

// ROC
TGraph *finePtBin_mass_curve[3][nFineBins];
TGraph *Lead_CA12_mass_curve[3][nPtBins];
TGraph *Lead_CA12_pt_curve[3];
TGraph *Lead_CA12_scaled_pt_curve[3];
TGraph *pTweights_curve;
TGraph *PtReweight_curve;


TCanvas * c1[3][nPtBins]; 
TCanvas * c3[3][nFineBins];
TCanvas * c2[nPtBins];
TPad *pad1[nPtBins];
TPad *pad2[nPtBins];


void defineStrings(/*TString *AlgoList, TString *binLabel,*/ TString *pTbins, TString *finePtBins){//, std::map<std::string, std::string> subjetMap){

  pTbins[0]="inclusive";
  pTbins[1]="gtr1000";
  pTbins[2]="500to1000";
  pTbins[3]="1000to1500";
  pTbins[4]="1500to2000";
  pTbins[5]="gtr2000";
  ptrange.push_back(std::make_pair(0,10000));
  ptrange.push_back(std::make_pair(1000,10000));
  ptrange.push_back(std::make_pair(500,1000));
  ptrange.push_back(std::make_pair(1000,1500));
  ptrange.push_back(std::make_pair(1500,2000));
  ptrange.push_back(std::make_pair(2000,10000));
  
  finePtBins[0]="0to250";
  finePtBins[1]="250to500";
  finePtBins[2]="500to750";
  finePtBins[3]="750to1000";
  finePtBins[4]="1000to1250";
  finePtBins[5]="1250to1500";
  finePtBins[6]="1500to1750";
  finePtBins[7]="1750to2000";
  finePtBins[8]="2000to2250";
  finePtBins[9]="2250to2500";
  finePtBins[10]="2500to2750";
  finePtBins[11]="2750to3000";

  pTbinsN[0]="Inclusive p_{T}^{CA12}";
  pTbinsN[1]="p_{T}^{CA12} > 1000 GeV";
  pTbinsN[2]="500 < p_{T}^{CA12} < 1000 GeV";
  pTbinsN[3]="1000 < p_{T}^{CA12} < 1500 GeV";
  pTbinsN[4]="1500 < p_{T}^{CA12} < 2000 GeV";
  pTbinsN[5]="p_{T}^{CA12} > 2000 GeV";
  
  pTFinebinsN[0]="0 < p_{T}^{CA12} < 250 GeV";
  pTFinebinsN[1]="250 < p_{T}^{CA12} < 500 GeV";
  pTFinebinsN[2]="500 < p_{T}^{CA12} < 750 GeV";
  pTFinebinsN[3]="750 < p_{T}^{CA12} < 1000 GeV";
  pTFinebinsN[4]="1000 < p_{T}^{CA12} < 1250 GeV";
  pTFinebinsN[5]="1250 < p_{T}^{CA12} < 1500 GeV";
  pTFinebinsN[6]="1500 < p_{T}^{CA12} < 1750 GeV";
  pTFinebinsN[7]="1750 < p_{T}^{CA12} < 2000 GeV";
  pTFinebinsN[8]="2000 < p_{T}^{CA12} < 2250 GeV";
  pTFinebinsN[9]="2250 < p_{T}^{CA12} < 2500 GeV";
  pTFinebinsN[10]="2500 < p_{T}^{CA12} < 2750 GeV";
  pTFinebinsN[11]="2750 < p_{T}^{CA12} < 3000 GeV";



}


static inline std::string &rtrim(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
  return s;
}


//store all the weights for different runNumbers:
std::map<long, float> k_factors;
std::map<long, float> filt_eff;
std::map<long, float> xs;
//std::map<int, float> ;




std::map<int, std::vector<float> *> var_E_vec;
std::map<int, std::vector<float> *> var_pt_vec;
std::map<int, std::vector<float> *> var_m_vec;
std::map<int, std::vector<float> *> var_eta_vec;
std::map<int, std::vector<float> *> var_phi_vec;
std::map<int, std::vector<float> *> var_emfrac_vec;
std::map<int, std::vector<float> *> var_Tau1_vec;
std::map<int, std::vector<float> *> var_Tau2_vec;
std::map<int, std::vector<float> *> var_Tau3_vec;
std::map<int, std::vector<float> *> var_WIDTH_vec;
std::map<int, std::vector<float> *> var_SPLIT12_vec;
std::map<int, std::vector<float> *> var_SPLIT23_vec;
std::map<int, std::vector<float> *> var_SPLIT34_vec;
std::map<int, std::vector<float> *> var_Dip12_vec;
std::map<int, std::vector<float> *> var_Dip13_vec;
std::map<int, std::vector<float> *> var_Dip23_vec;
std::map<int, std::vector<float> *> var_DipExcl12_vec;
std::map<int, std::vector<float> *> var_PlanarFlow_vec;
std::map<int, std::vector<float> *> var_Angularity_vec;
std::map<int, std::vector<float> *> var_QW_vec;
std::map<int, std::vector<float> *> var_PullMag_vec;
std::map<int, std::vector<float> *> var_PullPhi_vec;
std::map<int, std::vector<float> *> var_Pull_C00_vec;
std::map<int, std::vector<float> *> var_Pull_C01_vec;
std::map<int, std::vector<float> *> var_Pull_C10_vec;
std::map<int, std::vector<float> *> var_Pull_C11_vec;

std::map<int, std::vector<float> *> var_TauWTA1_vec; 
std::map<int, std::vector<float> *> var_TauWTA2_vec; 
std::map<int, std::vector<float> *> var_TauWTA3_vec; 
std::map<int, std::vector<float> *> var_TauWTA2TauWTA1_vec; 
std::map<int, std::vector<float> *> var_ZCUT12_vec;
std::map<int, std::vector<float> *> var_ZCUT23_vec;
std::map<int, std::vector<float> *> var_ZCUT34_vec;

std::map<int, std::vector<float> *> var_ActiveArea_vec;
std::map<int, std::vector<float> *> var_Aplanarity_vec;
std::map<int, std::vector<float> *> var_Sphericity_vec;
std::map<int, std::vector<float> *> var_ThrustMaj_vec;
std::map<int, std::vector<float> *> var_ThrustMin_vec;
std::map<int, std::vector<float> *> var_VoronoiArea_vec;





// jet clusters
std::vector<float> * var_cl_n_vec;
std::vector<float> * var_cl_pt_vec;
std::vector<float> * var_cl_eta_vec;
std::vector<float> * var_cl_phi_vec;


std::vector<float> * var_YFilt_vec;

Float_t var_massFraction_vec;
Float_t var_ktycut2_vec;



// electrons in
vector<TLorentzVector> * var_electrons_vec;
std::vector<TLorentzVector> electrons;

std::vector<float> * var_electronPt_vec;
std::vector<float> * var_electronEta_vec;
std::vector<float> * var_electronPhi_vec;
//std::vector<float> * var_electronT_vec;
std::vector<float> * var_el_ptcone20_vec;
std::vector<float> * var_el_etcone20_vec;

// muons in
std::vector<TLorentzVector> * var_muons_vec;
std::vector<TLorentzVector> muons;

std::vector<float> * var_muonPt_vec;
std::vector<float> * var_muonEta_vec;
std::vector<float> * var_muonPhi_vec;
std::vector<float> * var_mu_ptcone20_vec;
std::vector<float> * var_mu_etcone20_vec;
std::vector<float> * var_mu_charge_vec;


std::map<int, std::vector<std::vector < int>  > * > var_constit_index;
std::map<int, std::vector < int>  * > var_constit_n;

std::vector<std::vector <int> > * subjet_index;

std::vector<float> * var_subjets_E_vec;
std::vector<float> * var_subjets_pt_vec;
std::vector<float> * var_subjets_m_vec;
std::vector<float> * var_subjets_eta_vec;
std::vector<float> * var_subjets_phi_vec;
//std::vector<float> var_massdrop_vec;
//std::vector<float> var_yt_vec;
//std::map<int, std::vector<float> > var_Tau21_vec;


// variables that we write out to the outfile

// store one value each for truth, topo, groomed
std::vector<float> var_E;
std::vector<float> var_pt;
std::vector<float> var_m;
std::vector<float> var_eta;
std::vector<float> var_phi;
std::vector<float> var_emfrac;
std::vector<float> var_Tau1;
std::vector<float> var_Tau2;
std::vector<float> var_Tau3;
std::vector<float> var_WIDTH;
std::vector<float> var_SPLIT12;
std::vector<float> var_SPLIT23;
std::vector<float> var_SPLIT34;
std::vector<float> var_Dip12;
std::vector<float> var_Dip13;
std::vector<float> var_Dip23;
std::vector<float> var_DipExcl12;
std::vector<float> var_PlanarFlow;
std::vector<float> var_Angularity;
std::vector<float> var_QW;
std::vector<float> var_PullMag;
std::vector<float> var_PullPhi;
std::vector<float> var_Pull_C00;
std::vector<float> var_Pull_C01;
std::vector<float> var_Pull_C10;
std::vector<float> var_Pull_C11;
std::vector<float> var_Tau21;
Float_t var_YFilt;




// leptons out
std::vector<TLorentzVector> var_leptons;
std::vector<Float_t> var_ptcone20;
std::vector<Float_t> var_etcone20;
Float_t var_mllj;
Float_t var_mll;
Float_t var_ptll;

// store leading jet pt, before any truth matching or selections
Float_t var_leadingJetPt;

Int_t var_isElectronEvent = false;

// muons out
/*
std::vector<TLorentzVector> var_muons;
std::vector<float> var_mu_ptcone20;
std::vector<float> var_mu_etcone20;
*/
std::vector<float> var_charge;

// these variables are only stored for the subjets of the groomed jets, so we don't need a vector
Float_t var_subjets_E;
Float_t var_subjets_pt;
Float_t var_subjets_m;
Float_t var_subjets_eta;
Float_t var_subjets_phi;
Float_t var_massdrop;
Float_t var_yt;

// these variables are only in Lily's samples, so we need to add a flag to say we are running on those samples
std::vector<Float_t> var_TauWTA1; 
std::vector<Float_t> var_TauWTA2; 
std::vector<Float_t> var_TauWTA3; 
std::vector<Float_t> var_TauWTA2TauWTA1; 
std::vector<Float_t> var_ZCUT12;
std::vector<Float_t> var_ZCUT23;
std::vector<Float_t> var_ZCUT34;


std::vector<Float_t> var_ActiveArea;
std::vector<Float_t> var_Aplanarity;
std::vector<Float_t> var_Sphericity;
std::vector<Float_t> var_ThrustMaj;
std::vector<Float_t> var_ThrustMin;
std::vector<Float_t> var_VoronoiArea;




std::vector<Float_t> var_QjetVol;
std::vector<Float_t> var_FoxWolfram20;
std::vector<Float_t> var_softdrop;
std::vector<Float_t> var_EEC_C1;
std::vector<Float_t> var_EEC_C2;
std::vector<Float_t> var_EEC_D1;
std::vector<Float_t> var_EEC_D2;
// store the weights for the samples
Float_t var_k_factor;
Float_t var_filter_eff;
Float_t var_xs;



bool subjetscalc;
bool subjetspre;


